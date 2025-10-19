import pandas as pd 
import Bio 
from bio import Entrez
from bio import SeqIO


Entrez.email = "arin.kathapurkar26@gmail.com"


df1 = pd.read_excel("Aptamers1.xlsx") #can eventually combine both datasets into 1 
df_cleanedcols = df[["Type of Nucleic Acid","Target", "Aptamer Sequence",]] #only keep important columns 
target_col = df["Target"]


keeplist = []


for target in target_col:
    with Entrez.esearch(db= "protein", term = target_name) as handle:
        h = Entrez.read(handle)
    ids = h.get("IDList", [])
    keepmask.append[len(ids) > 0]
    
df proteins = df_cleanedcols[keeplist].reset_index(drop = True)

for target in targetcol:
    #write something so that I can do this more efficiently --> dictionary that holds protein: List[ids]
    with Entrez.esearch(db= "protein", term = target_name) as handle:
        h = Entrez.read(handle)
    ids = h.get("IDList", [])
    top_id = ids[0]
    with fasta = Entrez.efetch(ds = "protein", id = top_id, rettype = "fasta", retmode = "text") as f:
        fastatext = f.read()




    
        
        
        
        
       
        
        
    
    
    
        
    
     