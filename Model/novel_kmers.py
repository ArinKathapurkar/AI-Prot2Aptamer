import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import r2_score, mean_absolute_error
import pickle
from collections import Counter, defaultdict
import random


df = pd.read_csv('/Users/akathapurkar/Coding/Prot2Aptamer/Datasets/CLEAN/FINAL_apt_Data.csv')
feature_map = pd.read_csv('/Users/akathapurkar/Coding/Prot2Aptamer/Datasets/CLEAN/generation_feature_mapping.csv')


input_features = feature_map[feature_map['is_input'] == True]['feature_name'].tolist()
output_features = feature_map[feature_map['is_target'] == True]['feature_name'].tolist()

X = df[input_features].fillna(0)
y = df[output_features].fillna(0)


X_train, X_test, y_train, y_test = train_test_split(
    X, y, test_size=0.2, random_state=42
)


model = RandomForestRegressor(
    n_estimators=100,      # Number of trees
    max_depth=15,          # How deep each tree can go
    min_samples_leaf=2,    # Minimum samples at leaf nodes
    random_state=42,       # For reproducibility
    n_jobs=-1              # Use all CPU cores
)


model.fit(X_train, y_train)

y_pred = model.predict(X_test)

r2_scores = []
for i in range(y_test.shape[1]):
    r2 = r2_score(y_test.iloc[:, i], y_pred[:, i])
    r2_scores.append(r2)
    print(f"{output_features[i]:30s} R² = {r2:.3f}")


avg_r2 = np.mean(r2_scores)
print(f"\nAverage R²: {avg_r2:.3f}")
print(f"{'✓ GOOD!' if avg_r2 > 0.3 else 'Needs improvement'}")


with open('my_aptamer_predictor.pkl', 'wb') as f:
    pickle.dump(model, f)

print("✓ Model saved!")



def extract_protein_features_from_sequence(protein_seq):
    seq = str(protein_seq).upper()
    length = len(seq)
    
    if length == 0:
        print("ERROR: Empty protein sequence!")
        return None
    
    aa_counts = Counter(seq)
    features = {'prot_length': length}
    
    # Common amino acids (must match your training features!)
    common_aa = ['L', 'A', 'G', 'V', 'E', 'S', 'I', 'K', 'R', 'D']
    for aa in common_aa:
        features[f'prot_{aa}_pct'] = aa_counts.get(aa, 0) / length
    
    # Property groups
    hydrophobic = ['A', 'V', 'I', 'L', 'M', 'F', 'W', 'P']
    polar = ['S', 'T', 'Y', 'N', 'Q', 'C']
    charged_pos = ['K', 'R', 'H']
    charged_neg = ['D', 'E']
    aromatic = ['F', 'W', 'Y']
    
    features['prot_hydrophobic_pct'] = sum(aa_counts.get(aa, 0) for aa in hydrophobic) / length
    features['prot_polar_pct'] = sum(aa_counts.get(aa, 0) for aa in polar) / length
    features['prot_charged_positive_pct'] = sum(aa_counts.get(aa, 0) for aa in charged_pos) / length
    features['prot_charged_negative_pct'] = sum(aa_counts.get(aa, 0) for aa in charged_neg) / length
    features['prot_aromatic_pct'] = sum(aa_counts.get(aa, 0) for aa in aromatic) / length
    
    # Hydropathy (Kyte-Doolittle scale)
    HYDROPATHY = {
        'A': 1.8, 'R': -4.5, 'N': -3.5, 'D': -3.5, 'C': 2.5,
        'Q': -3.5, 'E': -3.5, 'G': -0.4, 'H': -3.2, 'I': 4.5,
        'L': 3.8, 'K': -3.9, 'M': 1.9, 'F': 2.8, 'P': -1.6,
        'S': -0.8, 'T': -0.7, 'W': -0.9, 'Y': -1.3, 'V': 4.2
    }
    hydropathy_scores = [HYDROPATHY.get(aa, 0) for aa in seq]
    features['prot_hydropathy_mean'] = np.mean(hydropathy_scores)
    
    # Charge
    positive = aa_counts.get('K', 0) + aa_counts.get('R', 0) + aa_counts.get('H', 0)
    negative = aa_counts.get('D', 0) + aa_counts.get('E', 0)
    features['prot_net_charge'] = (positive - negative) / length
    
    # Secondary structure propensity
    helix_formers = ['A', 'E', 'L', 'M', 'Q', 'K', 'R', 'H']
    sheet_formers = ['V', 'I', 'Y', 'F', 'W', 'T']
    features['prot_helix_propensity'] = sum(aa_counts.get(aa, 0) for aa in helix_formers) / length
    features['prot_sheet_propensity'] = sum(aa_counts.get(aa, 0) for aa in sheet_formers) / length
    
    # Entropy
    entropy = -sum((count/length) * np.log2(count/length) 
                   for count in aa_counts.values() if count > 0)
    features['prot_entropy'] = entropy
    
    return features


def build_kmer_library(sequences, k=3):
    """Extract all k-mers from real DNA aptamer sequences"""
    kmer_counts = Counter()
    for seq in sequences:
        seq = str(seq).upper()
        for i in range(len(seq) - k + 1):
            kmer = seq[i:i+k]
            if all(base in 'ACGT' for base in kmer):
                kmer_counts[kmer] += 1
    total = sum(kmer_counts.values())
    kmer_probs = {kmer: count/total for kmer, count in kmer_counts.items()}
    return kmer_probs


def build_transition_matrix(sequences, k=3):
    """Build a k-mer transition matrix"""
    transitions = defaultdict(Counter)
    for seq in sequences:
        seq = str(seq).upper()
        for i in range(len(seq) - k):
            current_kmer = seq[i:i+k]
            next_kmer = seq[i+1:i+k+1]
            if all(base in 'ACGT' for base in current_kmer + next_kmer):
                transitions[current_kmer][next_kmer] += 1
    transition_probs = {}
    for kmer, next_kmers in transitions.items():
        total = sum(next_kmers.values())
        transition_probs[kmer] = {
            next_kmer: count/total 
            for next_kmer, count in next_kmers.items()
        }
    return transition_probs


def sample_nucleotide(a_pct, c_pct, g_pct, t_pct):
    """Sample a single DNA nucleotide based on target percentages"""
    total = a_pct + c_pct + g_pct + t_pct
    if total == 0:
        return random.choice(['A', 'C', 'G', 'T'])
    probs = np.array([a_pct, c_pct, g_pct, t_pct]) / total
    return np.random.choice(['A', 'C', 'G', 'T'], p=probs)


def generate_dna_sequence(predicted_features, kmer_library, transition_matrix,
                         feature_names, k=3):
    """Generate a DNA aptamer sequence that matches predicted features"""
    feature_dict = {name: val for name, val in zip(feature_names, predicted_features)}
    
    target_length = int(feature_dict['apt_length'])
    target_gc = feature_dict['apt_gc_content']
    target_a_pct = feature_dict['apt_A_pct']
    target_c_pct = feature_dict['apt_C_pct']
    target_g_pct = feature_dict['apt_G_pct']
    target_t_pct = feature_dict['apt_T_pct']
    
    target_length = max(10, min(target_length, 150))
    
    gc_tolerance = 0.15
    suitable_kmers = []
    for kmer, prob in kmer_library.items():
        kmer_gc = (kmer.count('G') + kmer.count('C')) / len(kmer)
        if abs(kmer_gc - target_gc) < gc_tolerance:
            suitable_kmers.append(kmer)
    
    if not suitable_kmers:
        suitable_kmers = list(kmer_library.keys())
    
    sequence = random.choice(suitable_kmers)
    
    max_attempts = target_length * 3
    attempts = 0
    
    while len(sequence) < target_length and attempts < max_attempts:
        attempts += 1
        current_kmer = sequence[-k:]
        
        if current_kmer in transition_matrix:
            next_kmers = transition_matrix[current_kmer]
            
            filtered_next = {}
            for next_kmer, prob in next_kmers.items():
                kmer_gc = (next_kmer.count('G') + next_kmer.count('C')) / len(next_kmer)
                if abs(kmer_gc - target_gc) < gc_tolerance:
                    filtered_next[next_kmer] = prob
            
            if filtered_next:
                next_kmers_list = list(filtered_next.keys())
                next_probs = np.array(list(filtered_next.values()))
                next_probs = next_probs / next_probs.sum()
                chosen_kmer = np.random.choice(next_kmers_list, p=next_probs)
                sequence += chosen_kmer[-1]
            else:
                sequence += sample_nucleotide(target_a_pct, target_c_pct, 
                                              target_g_pct, target_t_pct)
        else:
            sequence += sample_nucleotide(target_a_pct, target_c_pct, 
                                          target_g_pct, target_t_pct)
    
    sequence = sequence[:target_length]
    return sequence


def validate_sequence(sequence, predicted_features, feature_names):
    """Check if generated DNA sequence matches predicted features"""
    length = len(sequence)
    a_count = sequence.count('A')
    c_count = sequence.count('C')
    g_count = sequence.count('G')
    t_count = sequence.count('T')
    
    actual_features = {
        'length': length,
        'gc_content': (g_count + c_count) / length if length > 0 else 0,
        'a_pct': a_count / length if length > 0 else 0,
        'c_pct': c_count / length if length > 0 else 0,
        'g_pct': g_count / length if length > 0 else 0,
        't_pct': t_count / length if length > 0 else 0,
    }
    
    feature_dict = {name: val for name, val in zip(feature_names, predicted_features)}
    
    matches = {
        'length_diff': abs(actual_features['length'] - feature_dict['apt_length']),
        'gc_diff': abs(actual_features['gc_content'] - feature_dict['apt_gc_content']),
        'a_pct_diff': abs(actual_features['a_pct'] - feature_dict['apt_A_pct']),
        'c_pct_diff': abs(actual_features['c_pct'] - feature_dict['apt_C_pct']),
        'g_pct_diff': abs(actual_features['g_pct'] - feature_dict['apt_G_pct']),
        't_pct_diff': abs(actual_features['t_pct'] - feature_dict['apt_T_pct']),
    }
    
    quality_score = (
        matches['length_diff'] / 100 +
        matches['gc_diff'] * 10 +
        matches['a_pct_diff'] * 5 +
        matches['c_pct_diff'] * 5 +
        matches['g_pct_diff'] * 5 +
        matches['t_pct_diff'] * 5
    )
    
    return {
        'actual_features': actual_features,
        'matches': matches,
        'quality_score': quality_score
    }


def generate_aptamer_candidates(protein_features, model, kmer_library,
                                transition_matrix, feature_names,
                                n_candidates=10, k=3):
    """Generate multiple DNA aptamer candidates and return best ones"""
    predicted_features = model.predict(protein_features.reshape(1, -1))[0]
    
    candidates = []
    
    for i in range(n_candidates):
        sequence = generate_dna_sequence(
            predicted_features,
            kmer_library,
            transition_matrix,
            feature_names,
            k=k
        )
        
        validation = validate_sequence(sequence, predicted_features, feature_names)
        
        candidates.append({
            'sequence': sequence,
            'quality_score': validation['quality_score'],
            'validation': validation
        })
    
    candidates.sort(key=lambda x: x['quality_score'])
    
    return candidates

def generate_aptamers_for_new_protein(protein_sequence, protein_name="Unknown Protein",
                                     model=None, kmer_library=None, 
                                     transition_matrix=None, feature_names=None,
                                     n_candidates=10):
    """
    Generate aptamers for a completely new protein sequence
    
    Args:
        protein_sequence: String of amino acids (e.g., "MKLPQR...")
        protein_name: Name/description of the protein
        model: Trained model
        kmer_library: K-mer library from real aptamers
        transition_matrix: K-mer transitions
        feature_names: List of output feature names
        n_candidates: How many aptamers to generate
    
    Returns:
        List of generated aptamer candidates
    """
    print("\n" + "="*80)
    print(f"GENERATING APTAMERS FOR NEW PROTEIN: {protein_name}")
    print("="*80)
    
    # Extract features from the new protein sequence
    print("\n1. Extracting protein features from sequence...")
    protein_features_dict = extract_protein_features_from_sequence(protein_sequence)
    
    if protein_features_dict is None:
        return None
    
    # Convert to array in correct order
    protein_features_array = np.array([protein_features_dict[feat] for feat in input_features])
    
    print(f"   ✓ Sequence length: {len(protein_sequence)} amino acids")
    print(f"   ✓ Extracted {len(protein_features_array)} features")
    
    # Show some key features
    print(f"\n   Key protein properties:")
    print(f"      Net charge: {protein_features_dict['prot_net_charge']:.3f}")
    print(f"      Hydropathy: {protein_features_dict['prot_hydropathy_mean']:.3f}")
    print(f"      Helix propensity: {protein_features_dict['prot_helix_propensity']:.2%}")
    
    # Generate aptamers
    print(f"\n2. Generating {n_candidates} aptamer candidates...")
    candidates = generate_aptamer_candidates(
        protein_features_array,
        model,
        kmer_library,
        transition_matrix,
        feature_names,
        n_candidates=n_candidates,
        k=3
    )
    
    # Display results
    print(f"\n🏆 TOP 5 GENERATED APTAMERS FOR {protein_name}:")
    print("="*80)
    for i, candidate in enumerate(candidates[:5], 1):
        seq = candidate['sequence']
        score = candidate['quality_score']
        val = candidate['validation']
        
        print(f"\n#{i} (Quality Score: {score:.4f})")
        print(f"   Sequence: {seq}")
        print(f"   Length: {val['actual_features']['length']} bp")
        print(f"   GC Content: {val['actual_features']['gc_content']:.2%}")
        print(f"   Composition: A={val['actual_features']['a_pct']:.1%}, "
              f"C={val['actual_features']['c_pct']:.1%}, "
              f"G={val['actual_features']['g_pct']:.1%}, "
              f"T={val['actual_features']['t_pct']:.1%}")
    
    print("\n" + "="*80)
    print("✅ GENERATION COMPLETE!")
    print("="*80)
    
    return candidates


# ==============================================================================
# BUILD K-MER LIBRARY (Do this once)
# ==============================================================================

print("\n" + "="*80)
print("BUILDING K-MER LIBRARY")
print("="*80)

real_sequences = df['aptamer_sequence'].tolist()
print(f"\nBuilding k-mer library from {len(real_sequences)} real DNA aptamers...")
kmer_library_3 = build_kmer_library(real_sequences, k=3)
transition_matrix_3 = build_transition_matrix(real_sequences, k=3)
print(f"✓ Built library with {len(kmer_library_3)} unique 3-mers")
print(f"✓ Built transition matrix with {len(transition_matrix_3)} transitions")



print("\n" + "="*80)
print("EXAMPLE 1: TESTING WITH PROTEIN FROM TEST SET")
print("="*80)

test_idx = 5
protein_features = X_test.iloc[test_idx].values
protein_name = df.iloc[X_test.index[test_idx]]['target_name']

print(f"\n🎯 Generating aptamers for: {protein_name[:60]}...")

candidates_test = generate_aptamer_candidates(
    protein_features,
    model,
    kmer_library_3,
    transition_matrix_3,
    output_features,
    n_candidates=5,
    k=3
)

for i, candidate in enumerate(candidates_test, 1):
    seq = candidate['sequence']
    score = candidate['quality_score']
    print(f"\n#{i}: {seq} (Score: {score:.4f})")


new_protein_sequence = input("Enter a Custom Protein")
new_protein_name = "My Custom Test Protein"

candidates_new = generate_aptamers_for_new_protein(
    protein_sequence=new_protein_sequence,
    protein_name=new_protein_name,
    model=model,
    kmer_library=kmer_library_3,
    transition_matrix=transition_matrix_3,
    feature_names=output_features,
    n_candidates=5
)