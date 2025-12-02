"""
This script:
1. Reads aptamer_interactions.csv
2. Fetches protein sequences from UniProt 
3. Creates a complete training dataset
"""

import pandas as pd
import requests
import time
from typing import Dict, Optional
import argparse
import sys

class UniProtSequenceFetcher:
    def __init__(self, verbose=True):
        """Initialize the UniProt sequence fetcher"""
        self.base_url = "https://rest.uniprot.org/uniprotkb"
        self.session = requests.Session()
        self.verbose = verbose
        self.cache = {}  # Cache to avoid refetching
    
    def fetch_sequence_by_id(self, uniprot_id: str) -> Optional[Dict]:
        """Fetch protein sequence from UniProt using the UniProt ID"""
        # Check cache first
        if uniprot_id in self.cache:
            return self.cache[uniprot_id]
        
        try:
            url = f"{self.base_url}/{uniprot_id}.json"
            response = self.session.get(url, timeout=15)
            response.raise_for_status()
            
            data = response.json()
            
            # Extract relevant information
            result = {
                'uniprot_id': data.get('primaryAccession', uniprot_id),
                'protein_name': '',
                'organism': '',
                'sequence': '',
                'length': 0,
                'gene_name': ''
            }
            
            # Get protein name
            if 'proteinDescription' in data:
                rec_name = data['proteinDescription'].get('recommendedName', {})
                if 'fullName' in rec_name:
                    result['protein_name'] = rec_name['fullName'].get('value', '')
            
            # Get organism
            if 'organism' in data:
                result['organism'] = data['organism'].get('scientificName', '')
            
            # Get sequence
            if 'sequence' in data:
                result['sequence'] = data['sequence'].get('value', '')
                result['length'] = data['sequence'].get('length', 0)
            
            # Get gene name
            if 'genes' in data and len(data['genes']) > 0:
                gene_names = data['genes'][0].get('geneName', {})
                if 'value' in gene_names:
                    result['gene_name'] = gene_names['value']
            
            # Cache the result
            self.cache[uniprot_id] = result
            return result
            
        except requests.exceptions.RequestException as e:
            if self.verbose:
                print(f"    ✗ Network error fetching {uniprot_id}: {e}")
            return None
        except Exception as e:
            if self.verbose:
                print(f"    ✗ Error fetching {uniprot_id}: {e}")
            return None
    
    def fetch_batch(self, uniprot_ids: list, delay=0.3) -> Dict[str, Dict]:
        """Fetch multiple protein sequences with rate limiting"""
        results = {}
        total = len(uniprot_ids)
        
        for i, uniprot_id in enumerate(uniprot_ids, 1):
            if self.verbose:
                print(f"[{i}/{total}] Fetching {uniprot_id}...", end=" ", flush=True)
            
            protein_data = self.fetch_sequence_by_id(uniprot_id)
            
            if protein_data and protein_data['sequence']:
                results[uniprot_id] = protein_data
                if self.verbose:
                    print(f"✓ {protein_data['length']} aa")
            else:
                if self.verbose:
                    print(f"✗ Failed")
            
            # Rate limiting to be respectful to UniProt
            if i < total:
                time.sleep(delay)
        
        return results
    
    def process_interactions_csv(self, 
                                 input_file='aptamer_interactions.csv',
                                 output_file='aptamer_protein_traianing_data.csv',
                                 protein_only=True,
                                 limit=None,
                                 min_aptamer_length=10,
                                 max_aptamer_length=200):
        """Process the aptamer_interactions.csv file and create training dataset"""
        
        print("="*80)
        print("  APTAMER-PROTEIN TRAINING DATA GENERATOR")
        print("="*80)
        
        # Load the interactions dataset
        print(f"\nStep 1: Loading {input_file}...")
        try:
            df = pd.read_csv(input_file)
            print(f"✓ Loaded {len(df)} interactions")
        except Exception as e:
            print(f"✗ Error loading file: {e}")
            return None
        
        # Filter for proteins with UniProt IDs
        print(f"\nStep 2: Filtering data...")
        if protein_only:
            df = df[df['ligand_type'] == 'Protein']
            print(f"✓ Filtered to {len(df)} protein interactions")
        
        df = df[df['target_uniprot'].notna() & (df['target_uniprot'] != '')]
        print(f"✓ {len(df)} interactions have valid UniProt IDs")
        
        # Filter by aptamer length
        df['aptamer_length_calc'] = df['aptamer_sequence'].str.len()
        df = df[
            (df['aptamer_length_calc'] >= min_aptamer_length) &
            (df['aptamer_length_calc'] <= max_aptamer_length)
        ]
        print(f"✓ {len(df)} interactions have aptamers between {min_aptamer_length}-{max_aptamer_length} bp")
        
        # Get unique protein IDs
        unique_proteins = df['target_uniprot'].unique()
        print(f"✓ {len(unique_proteins)} unique proteins to fetch")
        
        if limit:
            unique_proteins = unique_proteins[:limit]
            df = df[df['target_uniprot'].isin(unique_proteins)]
            print(f"⚠ Limited to {limit} proteins for testing")
        
        # Fetch sequences
        print(f"\nStep 3: Fetching protein sequences from UniProt...")
        print(f"(This may take a while - ~{len(unique_proteins) * 0.3 / 60:.1f} minutes)")
        
        protein_sequences = self.fetch_batch(unique_proteins.tolist())
        
        print(f"\n✓ Successfully fetched {len(protein_sequences)}/{len(unique_proteins)} protein sequences")
        
        if len(protein_sequences) == 0:
            print("✗ No protein sequences fetched. Cannot create dataset.")
            return None
        
        # Merge sequences with aptamer data
        print(f"\nStep 4: Merging aptamer-protein pairs...")
        
        df['protein_sequence'] = df['target_uniprot'].map(
            lambda x: protein_sequences.get(x, {}).get('sequence', '')
        )
        df['protein_length'] = df['target_uniprot'].map(
            lambda x: protein_sequences.get(x, {}).get('length', 0)
        )
        df['protein_name_uniprot'] = df['target_uniprot'].map(
            lambda x: protein_sequences.get(x, {}).get('protein_name', '')
        )
        df['protein_organism_uniprot'] = df['target_uniprot'].map(
            lambda x: protein_sequences.get(x, {}).get('organism', '')
        )
        df['gene_name'] = df['target_uniprot'].map(
            lambda x: protein_sequences.get(x, {}).get('gene_name', '')
        )
        
        # Filter to complete pairs
        df_complete = df[
            (df['aptamer_sequence'].notna()) & 
            (df['aptamer_sequence'] != '') &
            (df['protein_sequence'].notna()) & 
            (df['protein_sequence'] != '')
        ].copy()
        
        # Clean aptamer sequences (remove RNA if needed, keep DNA)
        df_complete['aptamer_sequence_clean'] = df_complete['aptamer_sequence'].str.upper()
        
        # Rename for clarity
        df_complete['aptamer_length'] = df_complete['aptamer_length_calc']
        
        # Select columns for output
        output_df = df_complete[[
            'aptamer_id',
            'target_id',
            'target_uniprot',
            'target_name',
            'protein_name_uniprot',
            'gene_name',
            'organism',
            'protein_organism_uniprot',
            'aptamer_sequence',
            'aptamer_length',
            'protein_sequence',
            'protein_length',
            'binding_conditions',
            'reference_pubmed_id',
            'ligand_type',
            'interaction_present'
        ]].copy()
        
        # Save to CSV
        print(f"\nStep 5: Saving to {output_file}...")
        output_df.to_csv(output_file, index=False)
        print(f"✓ Saved {len(output_df)} complete aptamer-protein pairs")
        
        # Print statistics
        self.print_statistics(output_df)
        
        return output_df
    
    def print_statistics(self, df):
        """Print dataset statistics"""
        print("\n" + "="*80)
        print("  DATASET STATISTICS")
        print("="*80)
        
        print(f"\nData Overview:")
        print(f"  Total aptamer-protein pairs: {len(df)}")
        print(f"  Unique proteins: {df['target_uniprot'].nunique()}")
        print(f"  Unique aptamers: {df['aptamer_id'].nunique()}")
        
        print(f"\nAptamer Statistics:")
        print(f"  Length - Min: {df['aptamer_length'].min()} bp")
        print(f"  Length - Max: {df['aptamer_length'].max()} bp")
        print(f"  Length - Mean: {df['aptamer_length'].mean():.1f} bp")
        print(f"  Length - Median: {df['aptamer_length'].median():.1f} bp")
        
        print(f"\nProtein Statistics:")
        print(f"  Length - Min: {df['protein_length'].min()} aa")
        print(f"  Length - Max: {df['protein_length'].max()} aa")
        print(f"  Length - Mean: {df['protein_length'].mean():.1f} aa")
        print(f"  Length - Median: {df['protein_length'].median():.1f} aa")
        
        print(f"\nTop 10 Organisms:")
        for org, count in df['organism'].value_counts().head(10).items():
            print(f"  {org}: {count}")
        
        print(f"\nTop 10 Proteins (by number of aptamers):")
        protein_counts = df.groupby(['target_uniprot', 'target_name']).size().sort_values(ascending=False).head(10)
        for (uniprot, name), count in protein_counts.items():
            name_short = name[:55] + "..." if len(name) > 55 else name
            print(f"  {uniprot} - {name_short}: {count} aptamers")
        
        print("\n" + "="*80)
        print("  SAMPLE PAIRS")
        print("="*80)
        for idx in range(min(3, len(df))):
            row = df.iloc[idx]
            print(f"\nPair {idx+1}:")
            print(f"  Target: {row['target_name'][:60]}")
            print(f"  UniProt: {row['target_uniprot']}")
            print(f"  Organism: {row['organism']}")
            print(f"  Protein ({row['protein_length']} aa): {row['protein_sequence'][:60]}...")
            print(f"  Aptamer ({row['aptamer_length']} bp): {row['aptamer_sequence'][:60]}...")
        
        print("\n" + "="*80)


def main():
    parser = argparse.ArgumentParser(
        description='Create aptamer-protein training dataset from interactions CSV'
    )
    parser.add_argument(
        '--input', '-i',
        default='datasets/aptamer_interactions.csv',
        help='Input CSV file (default: aptamer_interactions.csv)'
    )
    parser.add_argument(
        '--output', '-o',
        default='aptamer_protein_training_data.csv',
        help='Output CSV file (default: aptamer_protein_training_data.csv)'
    )
    parser.add_argument(
        '--limit', '-l',
        type=int,
        default=None,
        help='Limit number of proteins to fetch (for testing)'
    )
    parser.add_argument(
        '--protein-only',
        action='store_true',
        default=True,
        help='Only process protein targets (default: True)'
    )
    parser.add_argument(
        '--min-length',
        type=int,
        default=10,
        help='Minimum aptamer length in bp (default: 10)'
    )
    parser.add_argument(
        '--max-length',
        type=int,
        default=200,
        help='Maximum aptamer length in bp (default: 200)'
    )
    
    args = parser.parse_args()
    
    print("""
╔══════════════════════════════════════════════════════════════════════════════╗
║             APTAMER-PROTEIN TRAINING DATASET CREATOR                         ║
║                                                                              ║
║  This script fetches protein sequences from UniProt and pairs them          ║
║  with aptamer sequences to create a machine learning training dataset.      ║
╚══════════════════════════════════════════════════════════════════════════════╝
""")
    
    fetcher = UniProtSequenceFetcher(verbose=True)
    
    output_df = fetcher.process_interactions_csv(
        input_file=args.input,
        output_file=args.output,
        protein_only=args.protein_only,
        limit=args.limit,
        min_aptamer_length=args.min_length,
        max_aptamer_length=args.max_length
    )
    
    if output_df is not None and len(output_df) > 0:
        print(f"\n✅ SUCCESS! Training dataset created:")
        print(f"   📄 File: {args.output}")
        print(f"   📊 Pairs: {len(output_df)}")
        print(f"   🧬 Proteins: {output_df['target_uniprot'].nunique()}")
        print(f"\n🚀 Ready for machine learning!")
    else:
        print(f"\n❌ Failed to create dataset. Check errors above.")
        sys.exit(1)


if __name__ == "__main__":
    main()
