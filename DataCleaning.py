import pandas as pd 
import Bio 
from bio import Entrez
from bio import SeqIO
import time

Entrez.email = "arin.kathapurkar26@gmail.com"


df1 = pd.read_excel("Aptamers1.xlsx") 
df2 = pd.read_excel("Aptamers2.xlsx", skiprows = 1)
combined = pd.concat([df1, df2], ignore_index = True)
df_cleaned_cols = combined[["Type of Nucleic Acid","Target", "Aptamer Sequence","Sequence Length"]]
df_nodupes = df_cleanedcols.drop_duplicates(subset = ['Target'])
df_nodupes.to_excel("Cleaned_Aptamer_Data.xlsx", index=False)


#screening function to check if any of the targets have keywords that might mean they are a protein
def screen_prot(target_names):
    non_protein_keywords = ['cell', 'bacteria', 'virus', 'dna',
                            'rna', 'aptamer', 'oligonucleotide', 'peptide']
    # protein_keywords = ['protein', 'enzyme', 'receptor', 'antibody', 'kinase', 'factor', 'ligase', 'polymerase']
    for target in target_names: 
        if any(keyword in target.lower() for keyword in non_protein_keywords):
            target_names.remove(target)
    return target_names

target_col = df_nodupes["Target"]
screened_targets = screen_prot(target_col)


keeplist = []

for target in target_col:
    
    with Entrez.esearch(db= "protein", term = target_name) as handle:
        h = Entrez.read(handle)
    ids = h.get("IDList", [])
    keeplist.append[len(ids) > 0]
    
df proteins = df_cleanedcols[keeplist].reset_index(drop = True)





'''
for target in targetcol:
    with Entrez.esearch(db= "protein", term = target_name) as handle:
        h = Entrez.read(handle)
    ids = h.get("IDList", [])
    top_id = ids[0]
    with fasta = Entrez.efetch(ds = "protein", id = top_id, rettype = "fasta", retmode = "text") as f:
        fastatext = f.read()
'''



    
        
        
        
        
       
        
        
    
    
    
        
    
     