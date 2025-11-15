import ssl, certifi, urllib.request
import os
import pandas as pd 
import Bio 
from Bio import Entrez
from Bio import SeqIO
import time

Entrez.email = "arin.kathapurkar26@gmail.com"
Entrez.api_key = os.getenv("NCBI_API_KEY")

def normalize_header(df):
    df.columns = (
        df.columns.astype(str)
        .str.replace(r"[\u00A0\u200B\u200E\u200F]", "", regex=True)  # remove invisible spaces
        .str.replace(r"\s+", " ", regex=True)                       # collapse whitespace
        .str.strip()                                                # trim whitespace
    )
    return df



df1 = pd.read_excel("Datasets/Aptamers1.xlsx") 
df2 = pd.read_excel("Datasets/Aptamers2.xlsx", skiprows = 1)
df1 = normalize_header(df1)

combined = pd.concat([df1, df2], ignore_index = True)


df_cleanedcols = combined[["Type of Nucleic Acid","Target", "Aptamer Sequence","Sequence Length"]]
df_nodupes = df_cleanedcols.drop_duplicates(subset = ['Target'])
df_nodupes.to_excel("Cleaned_Aptamer_Data.xlsx", index=False)


#screening function to check if any of the targets have keywords that might mean they are a protein
def screen_prot(target_names):
    mask = []
    non_protein_keywords = ['cell', 'bacteria', 'virus', 'dna',
                            'rna', 'aptamer', 'oligonucleotide', 'peptide']
    # protein_keywords = ['protein', 'enzyme', 'receptor', 'antibody', 'kinase', 'factor', 'ligase', 'polymerase']
    for target in target_names: 
        if any(keyword in target.lower() for keyword in non_protein_keywords):
            mask.append(False)
        else:
            mask.append(True)
    return mask

target_col = df_nodupes["Target"].fillna("").astype(str).tolist()
mask = screen_prot(target_col)
filtered = df_nodupes[mask] #removes screened proteins


keeplist = []

for target in target_col:
    
    with Entrez.esearch(db= "protein", term = target) as handle:
        h = Entrez.read(handle)
    ids = h.get("IDList", [])
    keeplist.append[len(ids) > 0]
    
df_proteins = df_cleanedcols[keeplist].reset_index(drop = True)
print(keeplist)
    
        
        
        
        
       
        
        
    
    
    
        
    
     