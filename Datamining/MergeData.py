"""
Combine Aptamer-Protein Training Datasets

Merges aptamer_protein_training_data.csv (from interactions) and 
excel_aptamer_training_data.csv (from Excel files) into one master dataset.

USAGE:
    python combine_datasets.py
"""

import pandas as pd
import numpy as np

def combine_datasets(file1='aptamer_protein_training_data.csv',
                    file2='excel_aptamer_training_data.csv',
                    output='MASTER_aptamer_dataset.csv'):
    """Combine two aptamer datasets into one unified format"""
    
    print("="*80)
    print("  COMBINING APTAMER-PROTEIN DATASETS")
    print("="*80)
    
    # Load datasets
    print(f"\nStep 1: Loading datasets...")
    df1 = pd.read_csv(file1)
    df2 = pd.read_csv(file2)
    print(f"  ✓ Dataset 1 (interactions): {len(df1)} pairs")
    print(f"  ✓ Dataset 2 (Excel): {len(df2)} pairs")
    
    # Standardize column names
    print(f"\nStep 2: Standardizing columns...")
    
    # Dataset 1 already has good column names, just rename a few
    df1_clean = df1.copy()
    df1_clean = df1_clean.rename(columns={
        'target_uniprot': 'protein_uniprot_id',
        'protein_organism_uniprot': 'protein_organism',
        'reference_pubmed_id': 'pubmed_id'
    })
    
    # Dataset 2 needs some additions to match
    df2_clean = df2.copy()
    
    # Add missing columns with defaults
    df2_clean['aptamer_id'] = ['Excel_' + str(i) for i in range(len(df2))]
    df2_clean['target_id'] = df2_clean['protein_uniprot_id']
    df2_clean['organism'] = df2_clean['protein_organism']
    df2_clean['ligand_type'] = 'Protein'
    df2_clean['interaction_present'] = 1
    df2_clean['pubmed_id'] = df2_clean['pubmed_link'].fillna('')
    
    # Select common columns in same order
    common_columns = [
        'aptamer_id',
        'target_id',
        'protein_uniprot_id',
        'target_name',
        'protein_name_uniprot',
        'gene_name',
        'organism',
        'protein_organism',
        'aptamer_sequence',
        'aptamer_length',
        'protein_sequence',
        'protein_length',
        'binding_conditions',
        'pubmed_id',
        'ligand_type',
        'interaction_present'
    ]
    
    # Make sure all columns exist in both dataframes
    for col in common_columns:
        if col not in df1_clean.columns:
            df1_clean[col] = ''
        if col not in df2_clean.columns:
            df2_clean[col] = ''
    
    # Select only common columns
    df1_final = df1_clean[common_columns]
    df2_final = df2_clean[common_columns]
    
    print(f"  ✓ Standardized to {len(common_columns)} columns")
    
    # Combine datasets
    print(f"\nStep 3: Combining datasets...")
    df_combined = pd.concat([df1_final, df2_final], ignore_index=True)
    print(f"  ✓ Combined: {len(df_combined)} total pairs")
    
    # Remove exact duplicates
    print(f"\nStep 4: Removing duplicates...")
    initial_count = len(df_combined)
    
    # Remove exact sequence duplicates
    df_combined = df_combined.drop_duplicates(
        subset=['aptamer_sequence', 'protein_sequence'], 
        keep='first'
    )
    
    removed = initial_count - len(df_combined)
    print(f"  ✓ Removed {removed} duplicate pairs")
    print(f"  ✓ Final dataset: {len(df_combined)} unique pairs")
    
    # Validate sequences
    print(f"\nStep 5: Validating sequences...")
    
    def is_valid_dna(seq):
        if pd.isna(seq) or seq == '':
            return False
        return set(str(seq).upper()).issubset({'A', 'T', 'C', 'G', 'U', 'N'})
    
    def is_valid_protein(seq):
        if pd.isna(seq) or seq == '':
            return False
        return set(str(seq).upper()).issubset('ACDEFGHIKLMNPQRSTVWYX')
    
    df_combined['valid_aptamer'] = df_combined['aptamer_sequence'].apply(is_valid_dna)
    df_combined['valid_protein'] = df_combined['protein_sequence'].apply(is_valid_protein)
    
    invalid_aptamers = (~df_combined['valid_aptamer']).sum()
    invalid_proteins = (~df_combined['valid_protein']).sum()
    
    print(f"  ✓ Valid aptamer sequences: {df_combined['valid_aptamer'].sum()}/{len(df_combined)}")
    print(f"  ✓ Valid protein sequences: {df_combined['valid_protein'].sum()}/{len(df_combined)}")
    
    if invalid_aptamers > 0:
        print(f"  ⚠ {invalid_aptamers} invalid aptamer sequences (will be kept but flagged)")
    if invalid_proteins > 0:
        print(f"  ⚠ {invalid_proteins} invalid protein sequences (will be kept but flagged)")
    
    # Add metadata columns
    df_combined['dataset_source'] = ''
    df_combined.loc[df_combined['aptamer_id'].str.startswith('Apta_'), 'dataset_source'] = 'interactions'
    df_combined.loc[df_combined['aptamer_id'].str.startswith('Excel_'), 'dataset_source'] = 'excel'
    
    # Save
    print(f"\nStep 6: Saving combined dataset...")
    df_combined.to_csv(output, index=False)
    print(f"  ✓ Saved to {output}")
    
    # Print statistics
    print_statistics(df_combined)
    
    return df_combined


def print_statistics(df):
    """Print comprehensive dataset statistics"""
    print("\n" + "="*80)
    print("  MASTER DATASET STATISTICS")
    print("="*80)
    
    print(f"\nOverview:")
    print(f"  Total aptamer-protein pairs: {len(df)}")
    print(f"  Unique proteins (UniProt ID): {df['protein_uniprot_id'].nunique()}")
    print(f"  Unique aptamer sequences: {df['aptamer_sequence'].nunique()}")
    print(f"  Unique protein sequences: {df['protein_sequence'].nunique()}")
    
    print(f"\nDataset sources:")
    for source, count in df['dataset_source'].value_counts().items():
        print(f"  {source}: {count} pairs")
    
    print(f"\nAptamer statistics:")
    print(f"  Length - Min: {df['aptamer_length'].min():.0f} bp")
    print(f"  Length - Max: {df['aptamer_length'].max():.0f} bp")
    print(f"  Length - Mean: {df['aptamer_length'].mean():.1f} bp")
    print(f"  Length - Median: {df['aptamer_length'].median():.1f} bp")
    
    print(f"\nProtein statistics:")
    print(f"  Length - Min: {df['protein_length'].min():.0f} aa")
    print(f"  Length - Max: {df['protein_length'].max():.0f} aa")
    print(f"  Length - Mean: {df['protein_length'].mean():.1f} aa")
    print(f"  Length - Median: {df['protein_length'].median():.1f} aa")
    
    print(f"\nTop 10 organisms:")
    for org, count in df['organism'].value_counts().head(10).items():
        if pd.notna(org) and org != '':
            print(f"  {org}: {count}")
    
    print(f"\nTop 10 most-studied proteins:")
    protein_counts = df.groupby(['protein_uniprot_id', 'target_name']).size().sort_values(ascending=False).head(10)
    for (uniprot, name), count in protein_counts.items():
        if pd.notna(name) and name != '':
            name_short = name[:50] + "..." if len(name) > 50 else name
            print(f"  {uniprot} - {name_short}: {count} aptamers")
    
    # Validation summary
    valid_pairs = (df['valid_aptamer'] & df['valid_protein']).sum()
    print(f"\nData quality:")
    print(f"  Complete valid pairs: {valid_pairs}/{len(df)} ({valid_pairs/len(df)*100:.1f}%)")
    
    print("\n" + "="*80)
    print("  SAMPLE PAIRS FROM MASTER DATASET")
    print("="*80)
    
    # Show samples from each source
    for source in df['dataset_source'].unique():
        if pd.notna(source) and source != '':
            print(f"\nFrom {source} dataset:")
            df_source = df[df['dataset_source'] == source]
            for idx in range(min(2, len(df_source))):
                row = df_source.iloc[idx]
                print(f"  {idx+1}. {row['target_name'][:50]}")
                print(f"     Protein ({row['protein_length']:.0f} aa): {row['protein_sequence'][:50]}...")
                print(f"     Aptamer ({row['aptamer_length']:.0f} bp): {row['aptamer_sequence'][:50]}...")
    
    print("\n" + "="*80)
    print(f"\n✅ Master dataset ready for training!")
    print(f"   📊 {len(df)} pairs from {df['protein_uniprot_id'].nunique()} proteins")
    print(f"   📄 File: MASTER_aptamer_dataset.csv")
    print("="*80)


def main():
    import argparse
    
    parser = argparse.ArgumentParser(description='Combine aptamer-protein datasets')
    parser.add_argument('--file1', default='aptamer_protein_training_data.csv',
                       help='First dataset (interactions)')
    parser.add_argument('--file2', default='excel_aptamer_training_data.csv',
                       help='Second dataset (Excel)')
    parser.add_argument('--output', default='MASTER_aptamer_dataset.csv',
                       help='Output combined dataset')
    
    args = parser.parse_args()
    
    print("""
╔══════════════════════════════════════════════════════════════════════════════╗
║                 COMBINE APTAMER-PROTEIN DATASETS                             ║
║                                                                              ║
║  Merges multiple datasets into one unified training file                    ║
╚══════════════════════════════════════════════════════════════════════════════╝
""")
    
    try:
        df = combine_datasets(
            file1=args.file1,
            file2=args.file2,
            output=args.output
        )
        
        print(f"\n🎉 SUCCESS! Combined dataset created with {len(df)} pairs!")
        
    except Exception as e:
        print(f"\n❌ Error: {e}")
        import traceback
        traceback.print_exc()


if __name__ == "__main__":
    main()