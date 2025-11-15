from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
from selenium.webdriver.chrome.options import Options
from selenium.webdriver.common.action_chains import ActionChains
from bs4 import BeautifulSoup
import pandas as pd
import sqlite3
import time
import re
from typing import List, Dict

class AptamerScraper:
    def __init__(self, db_path='aptamer_database.db', headless=False):
        """Initialize scraper with Selenium WebDriver"""
        self.db_path = db_path
        self.base_url = "https://www.aptagen.com/apta-index/"
        
        # Setup Chrome options
        chrome_options = Options()
        if headless:
            chrome_options.add_argument('--headless')
        chrome_options.add_argument('--no-sandbox')
        chrome_options.add_argument('--disable-dev-shm-usage')
        chrome_options.add_argument('--window-size=1920,1080')
        
        self.driver = webdriver.Chrome(options=chrome_options)
        self.setup_database()
    
    def setup_database(self):
        """Create SQLite database and tables"""
        conn = sqlite3.connect(self.db_path)
        cursor = conn.cursor()
        
        # Check if table exists and what columns it has
        cursor.execute("SELECT name FROM sqlite_master WHERE type='table' AND name='aptamers'")
        table_exists = cursor.fetchone()
        
        if table_exists:
            # Get existing columns
            cursor.execute("PRAGMA table_info(aptamers)")
            existing_columns = [row[1] for row in cursor.fetchall()]
            
            # Expected columns
            expected_columns = {
                'id': 'INTEGER PRIMARY KEY AUTOINCREMENT',
                'target_name': 'TEXT',
                'aptamer_type': 'TEXT',
                'target_category': 'TEXT',
                'kd_value': 'TEXT',
                'sequence': 'TEXT',
                'sequence_length': 'INTEGER',
                'binding_conditions': 'TEXT',
                'binding_temp': 'TEXT',
                'doi': 'TEXT',
                'source_url': 'TEXT UNIQUE',
                'date_scraped': 'TIMESTAMP DEFAULT CURRENT_TIMESTAMP'
            }
            
            # Check if schema matches
            missing_columns = set(expected_columns.keys()) - set(existing_columns)
            if missing_columns or 'target_name' not in existing_columns:
                print(f"⚠ Database schema mismatch detected.")
                print(f"   Missing columns: {missing_columns}")
                print(f"   Existing columns: {existing_columns}")
                print(f"   Recreating table (this will delete existing data)...")
                # Drop and recreate table
                cursor.execute("DROP TABLE IF EXISTS aptamers")
                conn.commit()
        
        # Create table with correct schema
        cursor.execute('''
            CREATE TABLE IF NOT EXISTS aptamers (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                target_name TEXT,
                aptamer_type TEXT,
                target_category TEXT,
                kd_value TEXT,
                sequence TEXT,
                sequence_length INTEGER,
                binding_conditions TEXT,
                binding_temp TEXT,
                doi TEXT,
                source_url TEXT UNIQUE,
                date_scraped TIMESTAMP DEFAULT CURRENT_TIMESTAMP
            )
        ''')
        
        conn.commit()
        conn.close()
        print(f"✓ Database initialized at {self.db_path}")
    
    def apply_filters_on_index(self):
        """Apply DNA and Protein filters on the index page"""
        print("Loading Apta-Index and applying filters...")
        self.driver.get(self.base_url)
        
        # Wait for page to fully load
        wait = WebDriverWait(self.driver, 15)
        wait.until(EC.presence_of_element_located((By.TAG_NAME, "body")))
        time.sleep(3)  # Additional wait for dynamic content
        
        try:
            print("  Looking for filter controls...")
            
            # Debug: Print page structure to help identify filter elements
            self.debug_page_structure()
            
            # Try to apply both Protein and DNA filters
            protein_clicked = False
            dna_clicked = False
            
            # Helper function to click a filter element
            def click_filter_element(elem, filter_name):
                """Helper to click a filter element with multiple strategies"""
                try:
                    # Scroll element into view
                    self.driver.execute_script("arguments[0].scrollIntoView({block: 'center'});", elem)
                    time.sleep(0.5)
                    
                    # Check if element is visible and clickable
                    if not elem.is_displayed() or not elem.is_enabled():
                        return False
                    
                    # Try multiple click strategies
                    try:
                        # First try normal click
                        elem.click()
                    except:
                        # Try JavaScript click
                        try:
                            self.driver.execute_script("arguments[0].click();", elem)
                        except:
                            # Try mouse events
                            ActionChains(self.driver).move_to_element(elem).click().perform()
                    
                    time.sleep(1.5)  # Wait for filter to apply
                    print(f"  ✓ Clicked {filter_name} filter")
                    return True
                except Exception as e:
                    return False
            
            # Strategy 1: Try dropdown/select filters first (common in modern UIs)
            try:
                selects = self.driver.find_elements(By.TAG_NAME, "select")
                for select in selects:
                    try:
                        options = select.find_elements(By.TAG_NAME, "option")
                        for option in options:
                            option_text = option.text.strip()
                            if 'Protein' in option_text and not protein_clicked:
                                select.click()
                                time.sleep(0.5)
                                option.click()
                                time.sleep(2)
                                print("  ✓ Selected Protein from dropdown")
                                protein_clicked = True
                            elif 'DNA' in option_text and not dna_clicked:
                                select.click()
                                time.sleep(0.5)
                                option.click()
                                time.sleep(2)
                                print("  ✓ Selected DNA from dropdown")
                                dna_clicked = True
                    except:
                        continue
            except:
                pass
            
            # Strategy 2: Look for filter buttons/checkboxes - try BOTH filters
            filter_configs = [
                {'name': 'Protein', 'clicked': protein_clicked, 'var': 'protein_clicked'},
                {'name': 'DNA', 'clicked': dna_clicked, 'var': 'dna_clicked'}
            ]
            
            # First, try to expand any collapsed filter sections
            try:
                expand_buttons = self.driver.find_elements(By.XPATH, 
                    "//*[contains(@class, 'collapse') or contains(@class, 'accordion') or contains(@class, 'filter-panel')]//button")
                for btn in expand_buttons[:3]:  # Try first 3 expand buttons
                    try:
                        if btn.is_displayed():
                            btn.click()
                            time.sleep(1)
                    except:
                        pass
            except:
                pass
            
            for filter_config in filter_configs:
                filter_name = filter_config['name']
                if filter_config['clicked']:
                    continue  # Already clicked
                
                filter_selectors = [
                    f"//button[contains(., '{filter_name}')]",
                    f"//label[contains(., '{filter_name}')]",
                    f"//input[@type='checkbox' and (following-sibling::text()[contains(., '{filter_name}')] or preceding-sibling::text()[contains(., '{filter_name}')])]",
                    f"//span[contains(., '{filter_name}') and (@role='button' or @class[contains(., 'filter') or contains(., 'btn')])]",
                    f"//div[contains(@class, 'filter')]//*[contains(text(), '{filter_name}')]",
                    f"//a[contains(., '{filter_name}') and contains(@class, 'filter')]",
                    # Case-insensitive
                    f"//*[contains(translate(text(), 'ABCDEFGHIJKLMNOPQRSTUVWXYZ', 'abcdefghijklmnopqrstuvwxyz'), '{filter_name.lower()}') and (self::button or self::label or self::a or self::span)]",
                    # Look for filter containers
                    f"//*[contains(@class, 'filter') or contains(@class, 'facet')]//*[contains(text(), '{filter_name}')]",
                    # Look for clickable divs or li elements
                    f"//div[contains(., '{filter_name}') and (@onclick or @class[contains(., 'clickable') or contains(., 'filter')])]",
                    f"//li[contains(., '{filter_name}') and (@onclick or @class[contains(., 'clickable') or contains(., 'filter')])]",
                ]
                
                for selector in filter_selectors:
                    try:
                        elements = self.driver.find_elements(By.XPATH, selector)
                        if elements:
                            print(f"  Found {len(elements)} potential {filter_name} filter elements")
                            for elem in elements:
                                # Skip if this element is already selected/active (check for active classes)
                                try:
                                    classes = elem.get_attribute('class') or ''
                                    parent_classes = elem.find_element(By.XPATH, './..').get_attribute('class') or ''
                                    if 'active' in classes.lower() or 'selected' in classes.lower() or 'active' in parent_classes.lower():
                                        print(f"  {filter_name} filter appears already active, skipping click")
                                        if filter_name == 'Protein':
                                            protein_clicked = True
                                        else:
                                            dna_clicked = True
                                        break
                                except:
                                    pass
                                
                                if click_filter_element(elem, filter_name):
                                    if filter_name == 'Protein':
                                        protein_clicked = True
                                    else:
                                        dna_clicked = True
                                    break
                            if (filter_name == 'Protein' and protein_clicked) or (filter_name == 'DNA' and dna_clicked):
                                break
                    except:
                        continue
            
            # Wait for filters to apply and content to load
            if protein_clicked or dna_clicked:
                print("  Waiting for filtered content to load...")
                time.sleep(5)  # Wait for AJAX/filtered results to load
                
                # Wait for page to be ready
                try:
                    wait.until(lambda driver: driver.execute_script("return document.readyState") == "complete")
                except:
                    pass
                
                time.sleep(2)  # Additional wait for dynamic content
                
                # Print status
                print(f"  Filter Status: Protein={protein_clicked}, DNA={dna_clicked}")
            else:
                print("  ⚠ Could not find or click filter controls")
                print("  Will try to filter during scraping instead")
            
        except Exception as e:
            print(f"  ⚠ Error applying filters: {e}")
            print("  Will filter during scraping instead")
    
    def debug_page_structure(self):
        """Debug helper to inspect page structure"""
        print("\n  === PAGE DEBUG INFO ===")
        print(f"  Page Title: {self.driver.title}")
        print(f"  Current URL: {self.driver.current_url}")
        
        # Find all elements with "Protein" or "DNA" text
        protein_elements = self.driver.find_elements(By.XPATH, "//*[contains(text(), 'Protein')]")
        dna_elements = self.driver.find_elements(By.XPATH, "//*[contains(text(), 'DNA')]")
        
        print(f"  Elements containing 'Protein': {len(protein_elements)}")
        print(f"  Elements containing 'DNA': {len(dna_elements)}")
        
        # Print first few protein elements for debugging
        if protein_elements:
            print("\n  Sample Protein elements:")
            for i, elem in enumerate(protein_elements[:5]):
                try:
                    tag = elem.tag_name
                    text = elem.text[:50] if elem.text else "No text"
                    classes = elem.get_attribute('class') or "No class"
                    print(f"    {i+1}. <{tag}> class='{classes}' text='{text}'")
                except:
                    pass
        
        # Count all links
        all_links = self.driver.find_elements(By.TAG_NAME, "a")
        aptamer_links = [l for l in all_links if l.get_attribute('href') and '/aptamers/' in l.get_attribute('href')]
        print(f"\n  Total links: {len(all_links)}")
        print(f"  Aptamer links: {len(aptamer_links)}")
        print("  === END DEBUG ===\n")
    
    def find_entry_urls(self, max_scroll=15):
        """Find all aptamer detail page URLs from the index"""
        # Wait a bit for any dynamic content to load
        time.sleep(2)
        
        # Scroll to load more entries (lazy loading)
        print("  Scrolling to load all entries...")
        last_height = self.driver.execute_script("return document.body.scrollHeight")
        scroll_count = 0
        
        for i in range(max_scroll):
            # Scroll down
            self.driver.execute_script("window.scrollTo(0, document.body.scrollHeight);")
            time.sleep(2)  # Wait for lazy-loaded content
            
            # Also try scrolling a bit more incrementally
            self.driver.execute_script("window.scrollBy(0, 500);")
            time.sleep(1)
            
            new_height = self.driver.execute_script("return document.body.scrollHeight")
            if new_height == last_height:
                scroll_count += 1
                if scroll_count >= 2:  # If height hasn't changed for 2 scrolls, we're done
                    break
            else:
                scroll_count = 0
            last_height = new_height
            print(f"    Scroll {i+1}/{max_scroll} (height: {new_height})...")
        
        # Scroll back to top to ensure all links are in DOM
        self.driver.execute_script("window.scrollTo(0, 0);")
        time.sleep(1)
        
        # Find all aptamer links - try multiple strategies
        entry_urls = []
        
        # Strategy 1: Find all links with /aptamers/ in href
        print("  Searching for aptamer links...")
        links = self.driver.find_elements(By.TAG_NAME, "a")
        print(f"    Found {len(links)} total links on page")
        
        for link in links:
            try:
                href = link.get_attribute('href')
                if href and '/aptamers/' in href:
                    # Exclude non-aptamer pages and ensure it's a valid detail page
                    if not any(x in href for x in ['/apta-index/', '/technologies/', '#']):
                        # Make sure it's a full URL to an aptamer detail page
                        if href.startswith('http') and href not in entry_urls:
                            entry_urls.append(href)
            except:
                continue
        
        # Strategy 2: If no links found, try looking for links in specific containers
        if len(entry_urls) == 0:
            print("  Trying alternative link detection...")
            try:
                # Look for cards, tiles, or list items that might contain links
                containers = self.driver.find_elements(By.XPATH, 
                    "//*[contains(@class, 'card') or contains(@class, 'tile') or contains(@class, 'item')]//a")
                for container_link in containers:
                    try:
                        href = container_link.get_attribute('href')
                        if href and '/aptamers/' in href and href.startswith('http'):
                            if href not in entry_urls:
                                entry_urls.append(href)
                    except:
                        continue
            except:
                pass
        
        # Strategy 3: Try to find links by looking for aptamer-related text
        if len(entry_urls) == 0:
            print("  Trying text-based link detection...")
            try:
                # Look for any clickable elements near aptamer-related content
                aptamer_elements = self.driver.find_elements(By.XPATH,
                    "//a[contains(@href, 'aptamer') or contains(., 'aptamer')]")
                for elem in aptamer_elements:
                    try:
                        href = elem.get_attribute('href')
                        if href and href.startswith('http') and '/aptamers/' in href:
                            if href not in entry_urls:
                                entry_urls.append(href)
                    except:
                        continue
            except:
                pass
        
        # Debug: Print page title and URL to verify we're on the right page
        print(f"  Current page: {self.driver.title}")
        print(f"  Current URL: {self.driver.current_url}")
        
        print(f"✓ Found {len(entry_urls)} aptamer entry URLs")
        if len(entry_urls) > 0:
            print(f"  Sample URLs: {entry_urls[:3]}")
        
        return entry_urls
    
    def parse_nucleotide_sequence(self, seq_text):
        """Convert nucleotide notation to standard DNA sequence"""
        # Remove modifications like SH-(CH2)6-, biotin-, etc.
        seq_text = re.sub(r'^[^d]+?-(?=d[ATGC]p)', '', seq_text)
        
        # Map nucleotide codes to bases
        nucleotide_map = {
            'dAp': 'A', 'dA': 'A',
            'dTp': 'T', 'dT': 'T',
            'dGp': 'G', 'dG': 'G',
            'dCp': 'C', 'dC': 'C',
            'rAp': 'A', 'rA': 'A',
            'rUp': 'U', 'rU': 'U',
            'rGp': 'G', 'rG': 'G',
            'rCp': 'C', 'rC': 'C'
        }
        
        sequence = []
        # Find all nucleotide codes
        for code in re.findall(r'd[ATGC]p?|r[AUGC]p?', seq_text):
            if code in nucleotide_map:
                sequence.append(nucleotide_map[code])
        
        if sequence:
            return ''.join(sequence)
        
        # Fallback: if sequence is already in standard format
        clean_seq = re.sub(r'[^ATGCU]', '', seq_text.upper())
        if len(clean_seq) >= 10:
            return clean_seq
        
        return None
    
    def extract_aptamer_details(self, url: str) -> Dict:
        """Extract all information from an aptamer detail page"""
        try:
            self.driver.get(url)
            time.sleep(2)
            
            soup = BeautifulSoup(self.driver.page_source, 'html.parser')
            page_text = soup.get_text()
            
            # Split into lines for structured parsing
            lines = [line.strip() for line in page_text.split('\n') if line.strip()]
            
            data = {
                'source_url': url,
                'target_name': None,
                'aptamer_type': None,
                'target_category': None,
                'kd_value': None,
                'sequence': None,
                'sequence_length': None,
                'binding_conditions': None,
                'binding_temp': None,
                'doi': None
            }
            
            # Try to extract from HTML structure first (more reliable)
            # Look for definition lists, tables, or div structures
            if not data['target_category']:
                category_keywords = ['category', 'target category', 'antigen/target category']
                for keyword in category_keywords:
                    # Look for any element containing the keyword
                    for elem in soup.find_all(['dt', 'th', 'label', 'strong', 'span', 'div', 'p', 'td']):
                        elem_text = elem.get_text().strip()
                        if keyword in elem_text.lower():
                            # Check if category is in the same element text
                            for cat in ['Protein', 'Cells', 'Tissue', 'Small Organic', 'Peptide', 'Other', 'Antibody']:
                                if cat in elem_text:
                                    # Extract just the category part
                                    data['target_category'] = cat
                                    break
                            
                            # If not in same element, check next sibling
                            if not data['target_category']:
                                next_sibling = elem.find_next_sibling(['dd', 'td', 'div', 'span', 'p'])
                                if next_sibling:
                                    sibling_text = next_sibling.get_text().strip()
                                    for cat in ['Protein', 'Cells', 'Tissue', 'Small Organic', 'Peptide', 'Other', 'Antibody']:
                                        if cat in sibling_text:
                                            data['target_category'] = cat
                                            break
                            
                            # Also check parent's siblings
                            if not data['target_category']:
                                parent = elem.parent
                                if parent:
                                    for sibling in parent.find_next_siblings():
                                        sibling_text = sibling.get_text().strip()
                                        for cat in ['Protein', 'Cells', 'Tissue', 'Small Organic', 'Peptide', 'Other', 'Antibody']:
                                            if cat == sibling_text or cat in sibling_text:
                                                data['target_category'] = cat
                                                break
                                        if data['target_category']:
                                            break
                            
                            # Check following elements in the document
                            if not data['target_category']:
                                for next_elem in elem.find_next_siblings():
                                    next_elem_text = next_elem.get_text().strip()
                                    for cat in ['Protein', 'Cells', 'Tissue', 'Small Organic', 'Peptide', 'Other', 'Antibody']:
                                        if cat in next_elem_text:
                                            data['target_category'] = cat
                                            break
                                    if data['target_category']:
                                        break
                            
                            if data['target_category']:
                                break
                    if data['target_category']:
                        break
            
            # Extract target category - look for "Antigen/Target Category:" or similar patterns
            # Pattern 1: Look for "Category:" or "Target Category:" followed by the category
            category_match = re.search(
                r'(?:Antigen/)?Target\s+Category\s*:?\s*(Protein|Cells|Tissue|Small\s+Organic|Peptide|Other|Antibody)',
                page_text, re.IGNORECASE
            )
            if category_match:
                data['target_category'] = category_match.group(1)
            
            # Pattern 2: Look for "Category:" label followed by value on same or next line
            if not data['target_category']:
                category_label_match = re.search(
                    r'Category\s*:?\s*(Protein|Cells|Tissue|Small\s+Organic|Peptide|Other|Antibody)',
                    page_text, re.IGNORECASE
                )
                if category_label_match:
                    data['target_category'] = category_label_match.group(1)
            
            # Pattern 3: Look in the lines for label-value pairs
            if not data['target_category']:
                for i, line in enumerate(lines):
                    # Look for lines containing "Category" and the category value
                    if re.search(r'category', line, re.IGNORECASE):
                        # Check if this line contains a category value
                        for cat in ['Protein', 'Cells', 'Tissue', 'Small Organic', 'Peptide', 'Other', 'Antibody']:
                            if cat in line:
                                data['target_category'] = cat
                                break
                        # If not, check next line
                        if not data['target_category'] and i + 1 < len(lines):
                            next_line = lines[i + 1]
                            for cat in ['Protein', 'Cells', 'Tissue', 'Small Organic', 'Peptide', 'Other', 'Antibody']:
                                if cat == next_line:
                                    data['target_category'] = cat
                                    break
                        if data['target_category']:
                            break
            
            # Extract target name - try HTML structure first (more reliable)
            if not data['target_name']:
                # Look for "Target:" in HTML elements
                for elem in soup.find_all(['dt', 'th', 'label', 'strong', 'span', 'div', 'p', 'td']):
                    elem_text = elem.get_text().strip()
                    if 'target' in elem_text.lower() and ':' in elem_text:
                        # Check if target name is in the same element
                        # Extract everything after "Target:" or "Target"
                        target_match = re.search(r'Target\s*:?\s*(.+)', elem_text, re.IGNORECASE)
                        if target_match:
                            target_value = target_match.group(1).strip()
                            # Remove "Target:" if it's still there
                            target_value = re.sub(r'^Target\s*:?\s*', '', target_value, flags=re.IGNORECASE).strip()
                            if target_value and target_value.lower() != 'target':
                                data['target_name'] = target_value
                                break
                        
                        # If not in same element, check next sibling
                        if not data['target_name']:
                            next_sibling = elem.find_next_sibling(['dd', 'td', 'div', 'span', 'p'])
                            if next_sibling:
                                sibling_text = next_sibling.get_text().strip()
                                if sibling_text and sibling_text.lower() != 'target':
                                    data['target_name'] = sibling_text
                                    break
                    
                    if data['target_name']:
                        break
            
            # Find the aptamer type (DNA, RNA, etc.)
            for i, line in enumerate(lines):
                if line in ['DNA', 'RNA', '2\'-F-RNA', 'Chimeric', 'L-DNA', 'L-RNA', 'Peptide']:
                    data['aptamer_type'] = line
                    # Next line should be target name (if not already found)
                    if not data['target_name'] and i + 1 < len(lines):
                        next_line = lines[i + 1]
                        # Skip if it's just "Target:" label
                        if next_line.lower() not in ['target', 'target:']:
                            data['target_name'] = next_line
                    break
            
            # If target name not found yet, look for it more broadly with regex
            if not data['target_name']:
                # Look for "Target:" label followed by value
                # Pattern 1: "Target: Protein Name" on same line
                target_match = re.search(
                    r'Target\s*:?\s+([^\n]+?)(?:\n|$|Antigen|Category|Kd|Affinity|Binding)',
                    page_text, re.IGNORECASE | re.MULTILINE
                )
                if target_match:
                    target_value = target_match.group(1).strip()
                    # Clean up - remove common trailing words
                    target_value = re.sub(r'\s*(?:Antigen|Category|Kd|Affinity|Binding).*$', '', target_value, flags=re.IGNORECASE).strip()
                    if target_value and target_value.lower() not in ['target', 'target:']:
                        data['target_name'] = target_value
                
                # Pattern 2: Look for line after "Target:" label
                if not data['target_name']:
                    for i, line in enumerate(lines):
                        if re.search(r'^Target\s*:?\s*$', line, re.IGNORECASE):
                            if i + 1 < len(lines):
                                next_line = lines[i + 1].strip()
                                if next_line and next_line.lower() not in ['target', 'target:']:
                                    data['target_name'] = next_line
                                    break
            
            # Extract Kd value - look for patterns like "Kd = 1.2 nM" or "Kd: 1.2 nM"
            if not data['kd_value']:
                kd_patterns = [
                    r'Kd\s*[=:]\s*([0-9.]+\s*[µunpfm]*M)',
                    r'Kd\s+([0-9.]+\s*[µunpfm]*M)',
                    r'([0-9.]+\s*[µunpfm]*M)\s*Kd',
                ]
                for pattern in kd_patterns:
                    kd_match = re.search(pattern, page_text, re.IGNORECASE)
                    if kd_match:
                        data['kd_value'] = kd_match.group(1).strip()
                    break
            
            # Extract binding conditions
            cond_match = re.search(r'(0\.\d+\s*M\s+Na[^.]+\.)', page_text)
            if not cond_match:
                cond_match = re.search(r'(\d+\s*mM\s+[^.]+buffer[^.]+\.)', page_text, re.IGNORECASE)
            if cond_match:
                data['binding_conditions'] = cond_match.group(1).strip()
            
            # Extract binding temperature
            temp_match = re.search(r'(\d+\s*°C)', page_text)
            if temp_match:
                data['binding_temp'] = temp_match.group(1)
            
            # Extract DOI
            doi_match = re.search(r'doi:(10\.\d+/[^\s]+)', page_text)
            if doi_match:
                data['doi'] = doi_match.group(1)
            
            # Extract sequence - look for the nucleotide notation
            seq_match = re.search(r"5['\"]?\s*([^3]+?)\s*['\"]?3", page_text)
            if seq_match:
                seq_text = seq_match.group(1)
                sequence = self.parse_nucleotide_sequence(seq_text)
                if sequence:
                    data['sequence'] = sequence
                    data['sequence_length'] = len(sequence)
            
            # Alternative: look for dXpdXp pattern directly
            if not data['sequence']:
                seq_pattern = re.search(r'(d[ATGC]p(?:d[ATGC]p)+)', page_text)
                if seq_pattern:
                    sequence = self.parse_nucleotide_sequence(seq_pattern.group(1))
                    if sequence:
                        data['sequence'] = sequence
                        data['sequence_length'] = len(sequence)
            
            # Debug: Print extracted values to help diagnose issues
            if not data['target_category']:
                print(f"      ⚠ Warning: Could not extract target_category")
                print(f"         Aptamer type: {data['aptamer_type']}")
                print(f"         Target name: {data['target_name']}")
            
            return data
            
        except Exception as e:
            print(f"    ✗ Error: {e}")
            return None
    
    def save_to_database(self, data: Dict, update_existing=True):
        """Save aptamer data to database
        
        Args:
            data: Dictionary containing aptamer data
            update_existing: If True, update existing records based on source_url
        """
        if not data:
            return False
        
        conn = sqlite3.connect(self.db_path)
        cursor = conn.cursor()
        
        try:
            if update_existing:
                # Check if record exists based on source_url
                cursor.execute('SELECT id FROM aptamers WHERE source_url = ?', (data['source_url'],))
                existing = cursor.fetchone()
                
                if existing:
                    # Update existing record
                    cursor.execute('''
                        UPDATE aptamers SET
                            target_name = ?,
                            aptamer_type = ?,
                            target_category = ?,
                            kd_value = ?,
                            sequence = ?,
                            sequence_length = ?,
                            binding_conditions = ?,
                            binding_temp = ?,
                            doi = ?,
                            date_scraped = CURRENT_TIMESTAMP
                        WHERE source_url = ?
                    ''', (
                        data['target_name'],
                        data['aptamer_type'],
                        data['target_category'],
                        data['kd_value'],
                        data['sequence'],
                        data['sequence_length'],
                        data['binding_conditions'],
                        data['binding_temp'],
                        data['doi'],
                        data['source_url']
                    ))
                    action = "Updated"
                else:
                    # Insert new record
                    cursor.execute('''
                        INSERT INTO aptamers (target_name, aptamer_type, target_category,
                                             kd_value, sequence, sequence_length,
                                             binding_conditions, binding_temp, doi, source_url)
                        VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                    ''', (
                        data['target_name'],
                        data['aptamer_type'],
                        data['target_category'],
                        data['kd_value'],
                        data['sequence'],
                        data['sequence_length'],
                        data['binding_conditions'],
                        data['binding_temp'],
                        data['doi'],
                        data['source_url']
                    ))
                    action = "Saved"
            else:
                # Try INSERT first, if it fails due to duplicate, return False
                cursor.execute('''
                    INSERT INTO aptamers (target_name, aptamer_type, target_category,
                                         kd_value, sequence, sequence_length,
                                         binding_conditions, binding_temp, doi, source_url)
                    VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                ''', (
                    data['target_name'],
                    data['aptamer_type'],
                    data['target_category'],
                    data['kd_value'],
                    data['sequence'],
                    data['sequence_length'],
                    data['binding_conditions'],
                    data['binding_temp'],
                    data['doi'],
                    data['source_url']
                ))
                action = "Saved"
            
            conn.commit()
            seq_info = f"Seq: {data['sequence_length']}bp" if data['sequence'] else "No sequence"
            print(f"    ✓ {action}: {data['target_name']} ({data['aptamer_type']}, {seq_info})")
            return True
        except sqlite3.IntegrityError:
            if not update_existing:
                print(f"    ⚠ Already in database (skipped)")
                return False
            else:
                # Should not happen with INSERT OR REPLACE, but handle it
                print(f"    ✗ Database error: Could not save/update")
                return False
        except sqlite3.OperationalError as e:
            error_msg = str(e)
            if 'no such column' in error_msg.lower() or 'target name' in error_msg.lower():
                print(f"    ✗ Database schema error: {error_msg}")
                print(f"    ⚠ Try deleting {self.db_path} and running again to recreate the database")
            else:
                print(f"    ✗ Database error: {error_msg}")
            return False
        except Exception as e:
            print(f"    ✗ Database error: {e}")
            import traceback
            traceback.print_exc()
            return False
        finally:
            conn.close()
    
    def scrape_all(self, limit=None, filter_dna=True, filter_protein=True):
        """Main scraping method"""
        print("="*70)
        print("  APTAGEN SCRAPER - DNA Aptamers for Protein Targets")
        print("="*70)
        
        # Apply filters and get URLs
        self.apply_filters_on_index()
        entry_urls = self.find_entry_urls()
        
        if limit:
            entry_urls = entry_urls[:limit]
            print(f"\n⚠ Testing mode: Limited to {limit} entries\n")
        
        print(f"Processing {len(entry_urls)} entries...\n")
        
        successful = 0
        failed = 0
        skipped = 0
        
        for i, url in enumerate(entry_urls, 1):
            print(f"[{i}/{len(entry_urls)}] {url.split('/')[-2][:40]}...")
            
            data = self.extract_aptamer_details(url)
            
            if data:
                # Apply filters
                skip = False
                skip_reasons = []
                
                # Check DNA filter
                if filter_dna:
                    aptamer_type = data.get('aptamer_type', '').strip() if data.get('aptamer_type') else ''
                    if not aptamer_type:
                        print(f"    ⊘ Skipped: No aptamer type found")
                        skip = True
                        skip_reasons.append("no aptamer type")
                    elif 'DNA' not in aptamer_type.upper():
                        print(f"    ⊘ Skipped: aptamer_type='{aptamer_type}', not DNA")
                        skip = True
                        skip_reasons.append(f"not DNA (found: {aptamer_type})")
                
                # Check Protein filter
                if filter_protein and not skip:  # Only check if not already skipped
                    target_category = data.get('target_category', '').strip() if data.get('target_category') else ''
                    if not target_category:
                        print(f"    ⊘ Skipped: No target category found")
                        skip = True
                        skip_reasons.append("no target category")
                    elif target_category != 'Protein':
                        print(f"    ⊘ Skipped: target_category='{target_category}', not Protein")
                        skip = True
                        skip_reasons.append(f"not Protein (found: {target_category})")
                
                if skip:
                    skipped += 1
                    if skip_reasons:
                        print(f"      Reasons: {', '.join(skip_reasons)}")
                elif self.save_to_database(data):
                    successful += 1
                else:
                    failed += 1
            else:
                print(f"    ✗ Failed to extract data")
                failed += 1
            
            time.sleep(1)  # Be respectful
        
        print(f"\n{'='*70}")
        print(f"  RESULTS")
        print(f"{'='*70}")
        print(f"  ✓ Successful: {successful}")
        print(f"  ✗ Failed: {failed}")
        print(f"  ⊘ Skipped (filtered): {skipped}")
        print(f"{'='*70}\n")
    
    def view_sequences(self, limit=10, target_name=None):
        """View aptamer sequences from the database"""
        conn = sqlite3.connect(self.db_path)
        
        if target_name:
            query = "SELECT target_name, sequence, sequence_length FROM aptamers WHERE sequence IS NOT NULL AND target_name LIKE ? LIMIT ?"
            df = pd.read_sql_query(query, conn, params=(f'%{target_name}%', limit))
        else:
            query = "SELECT target_name, sequence, sequence_length FROM aptamers WHERE sequence IS NOT NULL LIMIT ?"
            df = pd.read_sql_query(query, conn, params=(limit,))
        
        conn.close()
        
        if len(df) > 0:
            print(f"\n{'='*70}")
            print(f"APTAMER SEQUENCES ({len(df)} shown):")
            print(f"{'='*70}")
            for idx, row in df.iterrows():
                print(f"\n{row['target_name']} ({row['sequence_length']}bp):")
                print(f"  {row['sequence']}")
            return df
        else:
            print("⚠ No sequences found in database")
            return pd.DataFrame()
    
    def export_to_csv(self, output_file='aptamers.csv'):
        """Export database to CSV"""
        conn = sqlite3.connect(self.db_path)
        df = pd.read_sql_query("SELECT * FROM aptamers", conn)
        conn.close()
        
        if len(df) > 0:
            df.to_csv(output_file, index=False)
            print(f"✓ Data exported to {output_file}")
            print(f"✓ Total entries: {len(df)}")
            print(f"✓ With sequences: {df['sequence'].notna().sum()}")
            return df
        else:
            print("⚠ No data to export")
            return pd.DataFrame()
    
    def close(self):
        """Close the browser"""
        self.driver.quit()


# Usage
if __name__ == "__main__":
    scraper = AptamerScraper(headless=False)
    
    try:
        # Test with 30 entries to see if it works
        scraper.scrape_all(
            limit=30,
            filter_dna=True,
            filter_protein=True
        )
        
        # Export results
        df = scraper.export_to_csv()
        
        if len(df) > 0:
            print("\n" + "="*70)
            print("SAMPLE DATA:")
            print("="*70)
            # Show sample data including sequences
            sample_cols = ['target_name', 'aptamer_type', 'target_category', 'sequence_length', 'sequence']
            available_cols = [col for col in sample_cols if col in df.columns]
            print(df[available_cols].head(10))
            
            # Also print sequences separately for easier viewing
            sequences_df = df[df['sequence'].notna()][['target_name', 'sequence']]
            if len(sequences_df) > 0:
                print("\n" + "="*70)
                print("APTAMER SEQUENCES:")
                print("="*70)
                for idx, row in sequences_df.head(10).iterrows():
                    print(f"\n{row['target_name']}:")
                    print(f"  {row['sequence']}")
                if len(sequences_df) > 10:
                    print(f"\n... and {len(sequences_df) - 10} more sequences (see CSV file for all)")
        
    finally:
        scraper.close()