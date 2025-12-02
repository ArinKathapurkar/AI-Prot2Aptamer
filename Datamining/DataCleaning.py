"""
Protein Sequence Fetcher for Aptamer Database
Fetches protein sequences from UniProt and pairs them with aptamer sequences
"""

import sqlite3
import requests
import time
import re
from typing import Dict, List, Optional, Tuple
import pandas as pd

class ProteinSequenceFetcher:
    def __init__(self, db_path='aptamer_database.db'):
        """Initialize protein sequence fetcher"""
        self.db_path = db_path
        self.uniprot_api = "https://rest.uniprot.org/uniprotkb/search"
        self.session = requests.Session()
        self.setup_database()
    
    def setup_database(self):
        """Add protein sequence columns to database if they don't exist"""
        conn = sqlite3.connect(self.db_path)
        cursor = conn.cursor()
        
        # Check existing columns
        cursor.execute("PRAGMA table_info(aptamers)")
        existing_columns = [row[1] for row in cursor.fetchall()]
        
        # Add new columns if they don't exist
        new_columns = {
            'protein_sequence': 'TEXT',
            'protein_uniprot_id': 'TEXT',
            'protein_organism': 'TEXT',
            'protein_length': 'INTEGER',
            'protein_fetch_status': 'TEXT',
            'protein_fetch_date': 'TIMESTAMP'
        }
        
        for col_name, col_type in new_columns.items():
            if col_name not in existing_columns:
                try:
                    cursor.execute(f"ALTER TABLE aptamers ADD COLUMN {col_name} {col_type}")
                    print(f"  ✓ Added column: {col_name}")
                except sqlite3.OperationalError as e:
                    if 'duplicate column' not in str(e).lower():
                        print(f"  ⚠ Could not add {col_name}: {e}")
        
        conn.commit()
        conn.close()
        print(f"✓ Database schema updated\n")
    
    def export_targets_for_manual_lookup(self, output_file='proteins_to_lookup.csv'):
        """Export list of protein targets that need sequences"""
        conn = sqlite3.connect(self.db_path)
        
        query = '''
            SELECT id, target_name, target_category, aptamer_type
            FROM aptamers
            WHERE protein_sequence IS NULL OR protein_sequence = ''
            ORDER BY target_name
        '''
        
        df = pd.read_sql_query(query, conn)
        conn.close()
        
        if len(df) > 0:
            # Add empty columns for manual entry
            df['protein_uniprot_id'] = ''
            df['protein_sequence'] = ''
            df['protein_organism'] = ''
            df['notes'] = ''
            
            df.to_csv(output_file, index=False)
            print(f"✓ Exported {len(df)} proteins to {output_file}")
            print(f"  Fill in: protein_uniprot_id, protein_sequence, protein_organism")
            print(f"  Then use import_manual_sequences() to add them to the database")
            return df
        else:
            print("⚠ All proteins already have sequences")
            return pd.DataFrame()
    
    def import_manual_sequences(self, input_file='proteins_to_lookup.csv'):
        """Import manually filled protein sequences from CSV"""
        try:
            df = pd.read_csv(input_file)
            
            required_cols = ['id', 'protein_uniprot_id', 'protein_sequence', 'protein_organism']
            if not all(col in df.columns for col in required_cols):
                print(f"✗ CSV must contain columns: {required_cols}")
                return
            
            # Filter rows that have been filled in
            df_filled = df[
                (df['protein_sequence'].notna()) & 
                (df['protein_sequence'] != '') &
                (df['protein_uniprot_id'].notna()) &
                (df['protein_uniprot_id'] != '')
            ].copy()
            
            if len(df_filled) == 0:
                print("⚠ No filled protein sequences found in CSV")
                return
            
            conn = sqlite3.connect(self.db_path)
            cursor = conn.cursor()
            
            success_count = 0
            for _, row in df_filled.iterrows():
                try:
                    protein_length = len(row['protein_sequence'])
                    
                    cursor.execute('''
                        UPDATE aptamers 
                        SET protein_sequence = ?,
                            protein_uniprot_id = ?,
                            protein_organism = ?,
                            protein_length = ?,
                            protein_fetch_status = 'manual',
                            protein_fetch_date = CURRENT_TIMESTAMP
                        WHERE id = ?
                    ''', (
                        row['protein_sequence'],
                        row['protein_uniprot_id'],
                        row['protein_organism'],
                        protein_length,
                        row['id']
                    ))
                    success_count += 1
                    print(f"  ✓ Updated ID {row['id']}: {row.get('target_name', 'Unknown')} ({protein_length} aa)")
                except Exception as e:
                    print(f"  ✗ Error updating ID {row['id']}: {e}")
            
            conn.commit()
            conn.close()
            
            print(f"\n✓ Imported {success_count}/{len(df_filled)} protein sequences")
            
        except Exception as e:
            print(f"✗ Error importing: {e}")
    
    def export_training_data(self, output_file='aptamer_training_data.csv', min_aptamer_length=10):
        """Export paired aptamer-protein sequences for ML training"""
        conn = sqlite3.connect(self.db_path)
        
        query = f'''
            SELECT 
                id,
                target_name,
                aptamer_type,
                target_category,
                sequence as aptamer_sequence,
                sequence_length as aptamer_length,
                protein_sequence,
                protein_length,
                protein_uniprot_id,
                protein_organism,
                kd_value,
                binding_conditions,
                binding_temp,
                protein_fetch_status
            FROM aptamers
            WHERE 
                sequence IS NOT NULL 
                AND sequence != ''
                AND protein_sequence IS NOT NULL 
                AND protein_sequence != ''
                AND sequence_length >= {min_aptamer_length}
            ORDER BY id
        '''
        
        df = pd.read_sql_query(query, conn)
        conn.close()
        
        if len(df) > 0:
            df.to_csv(output_file, index=False)
            print(f"\n✓ Training data exported to {output_file}")
            print(f"✓ Total paired sequences: {len(df)}")
            print(f"✓ Unique proteins: {df['protein_uniprot_id'].nunique()}")
            
            organisms = df['protein_organism'].value_counts()
            if len(organisms) > 0:
                print(f"✓ Top organisms: {dict(organisms.head(5))}")
            
            print(f"\nData statistics:")
            print(f"  - Aptamer length: {df['aptamer_length'].min()}-{df['aptamer_length'].max()} bp (mean: {df['aptamer_length'].mean():.1f})")
            print(f"  - Protein length: {df['protein_length'].min()}-{df['protein_length'].max()} aa (mean: {df['protein_length'].mean():.1f})")
            
            # Show sample
            print(f"\nSample paired sequences:")
            for idx in range(min(3, len(df))):
                row = df.iloc[idx]
                print(f"\n{idx+1}. {row['target_name']}")
                print(f"   Protein: {row['protein_sequence'][:60]}... ({row['protein_length']} aa)")
                print(f"   Aptamer: {row['aptamer_sequence'][:60]}... ({row['aptamer_length']} bp)")
            
            return df
        else:
            print("⚠ No complete paired sequences found")
            return pd.DataFrame()
    
    def show_summary(self):
        """Show summary of database status"""
        conn = sqlite3.connect(self.db_path)
        cursor = conn.cursor()
        
        cursor.execute("SELECT COUNT(*) FROM aptamers")
        total = cursor.fetchone()[0]
        
        cursor.execute("SELECT COUNT(*) FROM aptamers WHERE protein_sequence IS NOT NULL AND protein_sequence != ''")
        with_protein = cursor.fetchone()[0]
        
        cursor.execute("SELECT protein_fetch_status, COUNT(*) FROM aptamers GROUP BY protein_fetch_status")
        status_counts = dict(cursor.fetchall())
        
        cursor.execute("""
            SELECT COUNT(*) FROM aptamers 
            WHERE sequence IS NOT NULL AND sequence != '' 
            AND protein_sequence IS NOT NULL AND protein_sequence != ''
        """)
        complete_pairs = cursor.fetchone()[0]
        
        conn.close()
        
        print("\n" + "="*70)
        print("  DATABASE SUMMARY")
        print("="*70)
        print(f"Total aptamers: {total}")
        print(f"With protein sequences: {with_protein} ({with_protein/total*100:.1f}%)" if total > 0 else "With protein sequences: 0")
        print(f"Complete aptamer-protein pairs: {complete_pairs} ({complete_pairs/total*100:.1f}%)" if total > 0 else "Complete pairs: 0")
        print(f"\nFetch status breakdown:")
        for status, count in status_counts.items():
            status_label = status if status else 'Not attempted'
            print(f"  - {status_label}: {count}")
        print("="*70)


if __name__ == "__main__":
    fetcher = ProteinSequenceFetcher(db_path='aptamer_database.db')
    
    print("""
╔══════════════════════════════════════════════════════════════════════╗
║         PROTEIN SEQUENCE FETCHER FOR APTAMER DATABASE                ║
╚══════════════════════════════════════════════════════════════════════╝

This tool helps you fetch protein sequences to pair with your aptamers.

WORKFLOW:
1. Export targets that need protein sequences
2. Manually lookup sequences on UniProt (https://www.uniprot.org)
3. Import the filled CSV back into the database
4. Export training data for your ML model

""")
    
    # Show current status
    fetcher.show_summary()
    
    # Export proteins that need sequences
    print("\nStep 1: Exporting proteins that need sequences...")
    df = fetcher.export_targets_for_manual_lookup('proteins_to_lookup.csv')
    
    if len(df) > 0:
        print("\n" + "="*70)
        print("NEXT STEPS:")
        print("="*70)
        print("1. Open 'proteins_to_lookup.csv'")
        print("2. For each protein, search on https://www.uniprot.org")
        print("3. Fill in these columns:")
        print("   - protein_uniprot_id: e.g., P00734")
        print("   - protein_sequence: Copy the full amino acid sequence")
        print("   - protein_organism: e.g., Homo sapiens")
        print("4. Save the CSV file")
        print("5. Run: fetcher.import_manual_sequences('proteins_to_lookup.csv')")
        print("6. Then: fetcher.export_training_data()")
        print("="*70)
        
        print("\nShowing first 5 proteins to lookup:")
        print(df[['id', 'target_name', 'target_category']].head())