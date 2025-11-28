import time

protseq = "ALAARG"

def splitn(text, n):
    result = []
    for i in range (0, len(text), n):
        result.append(text[i:i+n])
    return result

protein_arr = splitn(protseq, 3)

#print(protein_arr)

#list of 3-letter amino acid codes, think about one-letter stuff 
aa_dict = {
    "ALA": 1,  # Alanine
    "ARG": 2,  # Arginine
    "ASN": 3,  # Asparagine
    "ASP": 4,  # Aspartic Acid
    "CYS": 5,  # Cysteine
    "GLU": 6,  # Glutamic Acid
    "GLN": 7,  # Glutamine
    "GLY": 8,  # Glycine
    "HIS": 9,  # Histidine
    "ILE": 10, # Isoleucine
    "LEU": 11, # Leucine
    "LYS": 12, # Lysine
    "MET": 13, # Methionine
    "PHE": 14, # Phenylalanine
    "PRO": 15, # Proline
    "SER": 16, # Serine
    "THR": 17, # Threonine
    "TRP": 18, # Tryptophan
    "TYR": 19, # Tyrosine
    "VAL": 20  # Valine
}


#turn amino acid codes into numbers
def tokenize(list):
    tokenized = ""
    for i in list:
        if i in aa_dict:
            tokenized = tokenized + str(aa_dict[i])
        else: 
            return "aa does not exist"
    return tokenized


token_seq = tokenize(protein_arr)

print(token_seq)

   

    


    
    


