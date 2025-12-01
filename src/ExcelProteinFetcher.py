"""
Process Excel Aptamer Files Using Existing UniProtSequenceFetcher

This script reuses your existing UniProtSequenceFetcher class to:
1. Load Aptamers1.xlsx and Aptamers2.xlsx
2. Clean and filter the data
3. Fetch protein sequences using your existing code
4. Create training dataset

USAGE:
    python process_excel_with_uniprot.py
"""

import pandas as pd
import re
import sys
import argparse

# Import your existing fetcher class
from UniProtLookup import UniProtSequenceFetcher


def clean_aptamer_sequence(seq):
    """Clean aptamer sequence - remove 5'/3' markers, keep only ATCGU"""
    if not seq or pd.isna(seq):
        return ""
    
    seq = str(seq).upper()
    seq = seq.replace("5'", "").replace("3'", "")
    seq = ''.join(c for c in seq if c in 'ATCGU')
    
    return seq


def is_likely_protein_target(target):
    """Determine if target is likely a protein (not a cell, molecule, etc.)"""
    if not target or pd.isna(target):
        return False
    
    target = str(target).lower()
    
    # Exclude non-protein targets
    exclude = ['cell', 'bacteria', 'virus particle', 'oocyst', 'spore',
               'dna', 'rna', 'atp', 'gtp', 'small molecule', 'dye', 
               'blue', 'reactive', 'cibacron', 'drug', 'toxin',
               'lactose', 'hpix', 'anep']
    
    if any(keyword in target for keyword in exclude):
        return False
    
    # Include likely proteins
    include = ['protein', 'enzyme', 'kinase', 'ase', 'factor',
               'receptor', 'antibody', 'immunoglobulin', 'albumin',
               'thrombin', 'fibrin', 'globin', 'transferrin',
               'integrin', 'cytokine', 'chemokine', 'interferon',
               'elastase', 'nucleolin', 'tubulin', 'fibronectin',
               'streptavidin', 'mucin', 'leptin', 'antigen',
               'fetoprotein', 'lipocalin', 'lactoferrin', 'rituximab']
    
    if any(keyword in target for keyword in include):
        return True
    
    # Ends with protein-like suffixes
    if any(target.endswith(suffix) for suffix in ['ase', 'in', 'gen', 'tor']):
        return True
    
    return False


def clean_target_name(target):
    """Clean target name for UniProt search"""
    if not target or pd.isna(target):
        return ""
    
    target = str(target).strip()
    
    # Remove everything in parentheses
    target = re.sub(r'\s*\([^)]*\)', '', target)
    
    # Remove comma and everything after
    target = re.sub(r',.*$', '', target)
    
    # Take first part if slash-separated
    target = target.split('/')[0]
    
    # Remove common prefixes
    target = re.sub(r'^(Recombinant|Human|recombinant)\s+', '', target, flags=re.IGNORECASE)
    
    return target.strip()


def search_uniprot_by_name(fetcher, target_name):
    """
    Search UniProt by protein name (not ID)
    This extends the fetcher to search by name instead of just ID
    """
    if not target_name:
        return None
    
    try:
        # Use UniProt search API
        search_url = "https://rest.uniprot.org/uniprotkb/search"
        
        # Clean name for search
        clean_name = clean_target_name(target_name)
        if not clean_name:
            return None
        
        # Build query - prefer human proteins
        query = f'(protein_name:{clean_name})'
        if 'human' not in clean_name.lower() and 'virus' not in clean_name.lower():
            query += ' AND (organism_name:human)'
        
        params = {
            'query': query,
            'format': 'json',
            'size': 1,
            'fields': 'accession,protein_name,organism_name,sequence,length,gene_names'
        }
        
        response = fetcher.session.get(search_url, params=params, timeout=15)
        response.raise_for_status()
        data = response.json()
        
        if 'results' in data and len(data['results']) > 0:
            entry = data['results'][0]
            
            result = {
                'uniprot_id': entry.get('primaryAccession', ''),
                'protein_name': '',
                'organism': '',
                'sequence': '',
                'length': 0,
                'gene_name': ''
            }
            
            # Get protein name
            if 'proteinDescription' in entry:
                rec_name = entry['proteinDescription'].get('recommendedName', {})
                if 'fullName' in rec_name:
                    result['protein_name'] = rec_name['fullName'].get('value', '')
            
            # Get organism
            if 'organism' in entry:
                result['organism'] = entry['organism'].get('scientificName', '')
            
            # Get sequence
            if 'sequence' in entry:
                result['sequence'] = entry['sequence'].get('value', '')
                result['length'] = entry['sequence'].get('length', 0)
            
            # Get gene name
            if 'genes' in entry and len(entry['genes']) > 0:
                gene_name_obj = entry['genes'][0].get('geneName', {})
                if 'value' in gene_name_obj:
                    result['gene_name'] = gene_name_obj['value']
            
            return result
    
    except Exception:
        pass
    
    return None


def process_excel_files(file1='Aptamers1.xlsx', 
                        file2='Aptamers2.xlsx',
                        output_file='excel_aptamer_training_data.csv',
                        limit=None):
    """Process Excel files and create training dataset"""
    
    print("="*80)
    print("  EXCEL APTAMER DATASET PROCESSOR")
    print("="*80)
    
    # Step 1: Load files
    print(f"\nStep 1: Loading Excel files...")
    df1 = pd.read_excel(file1)
    df2 = pd.read_excel(file2)
    print(f"  ✓ Loaded {len(df1)} from {file1}")
    print(f"  ✓ Loaded {len(df2)} from {file2}")
    
    # Step 2: Standardize and combine
    print(f"\nStep 2: Combining datasets...")
    df1 = df1.rename(columns={'Target': 'target'})
    df2 = df2.rename(columns={'Target ': 'target'})
    
    df = pd.concat([df1, df2], ignore_index=True)
    print(f"  ✓ Combined: {len(df)} entries")
    
    df = df.drop_duplicates(subset=['Aptamer Sequence'], keep='first')
    print(f"  ✓ After deduplication: {len(df)}")
    
    # Step 3: Filter for DNA
    print(f"\nStep 3: Filtering for DNA aptamers...")
    df = df[df['Type of Nucleic Acid'] == 'ssDNA']
    print(f"  ✓ DNA aptamers: {len(df)}")
    
    # Step 4: Filter for proteins
    print(f"\nStep 4: Filtering for protein targets...")
    df['is_protein'] = df['target'].apply(is_likely_protein_target)
    df = df[df['is_protein']]
    print(f"  ✓ Protein targets: {len(df)}")
    
    # Step 5: Clean sequences
    print(f"\nStep 5: Cleaning aptamer sequences...")
    df['aptamer_sequence'] = df['Aptamer Sequence'].apply(clean_aptamer_sequence)
    df['aptamer_length'] = df['aptamer_sequence'].str.len()
    
    # Filter by length
    df = df[(df['aptamer_length'] >= 10) & (df['aptamer_length'] <= 200)]
    print(f"  ✓ After length filter (10-200 bp): {len(df)}")
    
    # Step 6: Get unique targets
    print(f"\nStep 6: Identifying unique protein targets...")
    unique_targets = df['target'].unique()
    print(f"  ✓ {len(unique_targets)} unique targets")
    
    if limit:
        unique_targets = unique_targets[:limit]
        df = df[df['target'].isin(unique_targets)]
        print(f"  ⚠ Limited to {limit} targets for testing")
    
    # Step 7: Fetch protein sequences using YOUR existing class
    print(f"\nStep 7: Fetching protein sequences from UniProt...")
    print(f"(Estimated time: ~{len(unique_targets) * 0.3 / 60:.1f} minutes)\n")
    
    fetcher = UniProtSequenceFetcher(verbose=True)
    protein_map = {}
    
    success = 0
    failed = 0
    
    for i, target in enumerate(unique_targets, 1):
        print(f"[{i}/{len(unique_targets)}] {target[:50]}...", end=" ", flush=True)
        
        # Search by name (not ID)
        protein_data = search_uniprot_by_name(fetcher, target)
        
        if protein_data and protein_data['sequence']:
            protein_map[target] = protein_data
            print(f"✓ {protein_data['uniprot_id']} ({protein_data['length']} aa)")
            success += 1
        else:
            print(f"✗ Not found")
            failed += 1
        
        import time
        time.sleep(0.3)  # Rate limiting
    
    print(f"\n  ✓ Successfully fetched {success}/{len(unique_targets)} proteins")
    
    # Step 8: Merge with aptamer data
    print(f"\nStep 8: Merging aptamer-protein pairs...")
    
    df['protein_uniprot_id'] = df['target'].map(lambda x: protein_map.get(x, {}).get('uniprot_id', ''))
    df['protein_name_uniprot'] = df['target'].map(lambda x: protein_map.get(x, {}).get('protein_name', ''))
    df['protein_sequence'] = df['target'].map(lambda x: protein_map.get(x, {}).get('sequence', ''))
    df['protein_organism'] = df['target'].map(lambda x: protein_map.get(x, {}).get('organism', ''))
    df['protein_length'] = df['target'].map(lambda x: protein_map.get(x, {}).get('length', 0))
    df['gene_name'] = df['target'].map(lambda x: protein_map.get(x, {}).get('gene_name', ''))
    
    # Filter to complete pairs
    df_complete = df[
        (df['aptamer_sequence'] != '') &
        (df['protein_sequence'] != '')
    ].copy()
    
    print(f"  ✓ Complete pairs: {len(df_complete)}")
    
    # Step 9: Create output dataframe
    output_df = df_complete[[
        'target',
        'protein_uniprot_id',
        'protein_name_uniprot',
        'protein_sequence',
        'protein_organism',
        'protein_length',
        'gene_name',
        'aptamer_sequence',
        'aptamer_length',
        'Kd (nM)',
        'Affinity',
        'Binding Buffer/Conditions',
        'Year of Paper',
        'Link to PubMed Entry'
    ]].copy()
    
    # Rename for consistency
    output_df.columns = [
        'target_name',
        'protein_uniprot_id',
        'protein_name_uniprot',
        'protein_sequence',
        'protein_organism',
        'protein_length',
        'gene_name',
        'aptamer_sequence',
        'aptamer_length',
        'kd_nm',
        'affinity',
        'binding_conditions',
        'year',
        'pubmed_link'
    ]
    
    # Save
    output_df.to_csv(output_file, index=False)
    print(f"\n✓ Saved to {output_file}")
    
    # Print statistics
    print_statistics(output_df)
    
    return output_df


def print_statistics(df):
    """Print dataset statistics"""
    print("\n" + "="*80)
    print("  FINAL STATISTICS")
    print("="*80)
    print(f"Total pairs: {len(df)}")
    print(f"Unique proteins: {df['protein_uniprot_id'].nunique()}")
    print(f"Unique aptamers: {df['aptamer_sequence'].nunique()}")
    
    print(f"\nAptamer length:")
    print(f"  Min: {df['aptamer_length'].min()} bp")
    print(f"  Max: {df['aptamer_length'].max()} bp")
    print(f"  Mean: {df['aptamer_length'].mean():.1f} bp")
    
    print(f"\nProtein length:")
    print(f"  Min: {df['protein_length'].min()} aa")
    print(f"  Max: {df['protein_length'].max()} aa")
    print(f"  Mean: {df['protein_length'].mean():.1f} aa")
    
    print(f"\nTop 10 proteins:")
    top_proteins = df.groupby(['protein_uniprot_id', 'target_name']).size().sort_values(ascending=False).head(10)
    for (uniprot, name), count in top_proteins.items():
        print(f"  {uniprot} ({name[:40]}): {count} aptamers")
    
    print(f"\nYears: {df['year'].min():.0f}-{df['year'].max():.0f}")
    
    # Sample
    print("\n" + "="*80)
    print("  SAMPLE DATA")
    print("="*80)
    for idx in range(min(3, len(df))):
        row = df.iloc[idx]
        print(f"\n{idx+1}. {row['target_name']}")
        print(f"   UniProt: {row['protein_uniprot_id']}")
        print(f"   Protein: {row['protein_sequence'][:50]}... ({row['protein_length']} aa)")
        print(f"   Aptamer: {row['aptamer_sequence'][:50]}... ({row['aptamer_length']} bp)")
    
    print("\n" + "="*80)


def main():
    parser = argparse.ArgumentParser(
        description='Process Excel aptamer files using existing UniProtSequenceFetcher'
    )
    parser.add_argument('--file1', default='~/Coding/Prot2Aptamer/Datasets/Aptamers1.xlsx', help='First Excel file')
    parser.add_argument('--file2', default='Aptamers2.xlsx', help='Second Excel file')
    parser.add_argument('--output', default='excel_aptamer_training_data.csv', help='Output CSV')
    parser.add_argument('--limit', type=int, default=None, help='Limit proteins (for testing)')
    
    args = parser.parse_args()
    
    print("""
╔══════════════════════════════════════════════════════════════════════════════╗
║        EXCEL APTAMER PROCESSOR (Using Your Existing UniProt Fetcher)         ║
╚══════════════════════════════════════════════════════════════════════════════╝

This uses your existing UniProtSequenceFetcher class from UniProtLookup.py
to fetch protein sequences for the Excel datasets.
""")
    
    try:
        df = process_excel_files(
            file1=args.file1,
            file2=args.file2,
            output_file=args.output,
            limit=args.limit
        )
        
        if df is not None and len(df) > 0:
            print(f"\n✅ SUCCESS! Created {len(df)} aptamer-protein pairs")
            print(f"   📄 File: {args.output}")
            print(f"   🧬 Proteins: {df['protein_uniprot_id'].nunique()}")
            print(f"\n🚀 Ready for ML training!")
        else:
            print(f"\n❌ No complete pairs created")
            sys.exit(1)
    
    except Exception as e:
        print(f"\n❌ Error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()