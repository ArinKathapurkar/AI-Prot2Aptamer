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


new_protein_features = X_test.iloc[0:1]
predicted_aptamer_features = model.predict(new_protein_features)[0]
print("\nPredicted aptamer should have:")
for i, feature_name in enumerate(output_features):
    print(f"  {feature_name:30s} = {predicted_aptamer_features[i]:.3f}")

def build_kmer_library(sequences, k=3):
    """
    Extract all k-mers from real DNA aptamer sequences
    
    Args:
        sequences: List of DNA aptamer sequences (A, C, G, T only)
        k: K-mer size (default 3)
    
    Returns:
        Dictionary of k-mer probabilities
    """
    kmer_counts = Counter()
    
    for seq in sequences:
        seq = str(seq).upper()  # DNA only, no U→T conversion needed
        for i in range(len(seq) - k + 1):
            kmer = seq[i:i+k]
            if all(base in 'ACGT' for base in kmer):
                kmer_counts[kmer] += 1
    
    # Convert to probabilities
    total = sum(kmer_counts.values())
    kmer_probs = {kmer: count/total for kmer, count in kmer_counts.items()}
    
    return kmer_probs


def build_transition_matrix(sequences, k=3):
    """
    Build a k-mer transition matrix (which k-mers follow which)
    """
    transitions = defaultdict(Counter)
    
    for seq in sequences:
        seq = str(seq).upper()
        for i in range(len(seq) - k):
            current_kmer = seq[i:i+k]
            next_kmer = seq[i+1:i+k+1]
            
            if all(base in 'ACGT' for base in current_kmer + next_kmer):
                transitions[current_kmer][next_kmer] += 1
    
    # Convert to probabilities
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


def generate_dna_sequence(predicted_features, 
                         kmer_library, 
                         transition_matrix,
                         feature_names,
                         k=3):
    """
    Generate a DNA aptamer sequence that matches predicted features
    """
    # Extract key predicted features
    feature_dict = {name: val for name, val in zip(feature_names, predicted_features)}
    
    target_length = int(feature_dict['apt_length'])
    target_gc = feature_dict['apt_gc_content']
    target_a_pct = feature_dict['apt_A_pct']
    target_c_pct = feature_dict['apt_C_pct']
    target_g_pct = feature_dict['apt_G_pct']
    target_t_pct = feature_dict['apt_T_pct']
    
    # Ensure reasonable length
    target_length = max(10, min(target_length, 150))
    
    # Filter k-mers by GC content
    gc_tolerance = 0.15
    suitable_kmers = []
    for kmer, prob in kmer_library.items():
        kmer_gc = (kmer.count('G') + kmer.count('C')) / len(kmer)
        if abs(kmer_gc - target_gc) < gc_tolerance:
            suitable_kmers.append(kmer)
    
    if not suitable_kmers:
        suitable_kmers = list(kmer_library.keys())
    
    # Start sequence
    sequence = random.choice(suitable_kmers)
    
    # Extend using transition matrix
    max_attempts = target_length * 3
    attempts = 0
    
    while len(sequence) < target_length and attempts < max_attempts:
        attempts += 1
        current_kmer = sequence[-k:]
        
        if current_kmer in transition_matrix:
            next_kmers = transition_matrix[current_kmer]
            
            # Filter by GC content
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
    
    # Trim to exact length
    sequence = sequence[:target_length]
    return sequence


def validate_sequence(sequence, predicted_features, feature_names):
    """
    Check if generated DNA sequence matches predicted features
    """
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


def generate_aptamer_candidates(protein_features, 
                                model,
                                kmer_library,
                                transition_matrix,
                                feature_names,
                                n_candidates=10,
                                k=3):
    """
    Generate multiple DNA aptamer candidates and return best ones
    """
    # Predict aptamer features
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
    
    # Sort by quality
    candidates.sort(key=lambda x: x['quality_score'])
    
    return candidates

#GENERATING SEQUENCES

print("\n" + "="*80)
print("GENERATING DNA APTAMER SEQUENCES")
print("="*80)

# Get real DNA aptamer sequences
real_sequences = df['aptamer_sequence'].tolist()

# Build k-mer library
print("\nBuilding k-mer library from real DNA aptamers...")
kmer_library_3 = build_kmer_library(real_sequences, k=3)
transition_matrix_3 = build_transition_matrix(real_sequences, k=3)
print(f"✓ Built library with {len(kmer_library_3)} unique 3-mers")
print(f"✓ Built transition matrix with {len(transition_matrix_3)} transitions")

# Show some common k-mers
print("\nMost common 3-mers:")
top_kmers = sorted(kmer_library_3.items(), key=lambda x: x[1], reverse=True)[:10]
for kmer, freq in top_kmers:
    print(f"  {kmer}: {freq*100:.2f}%")

# Pick a test protein
test_idx = 5
protein_features = X_test.iloc[test_idx].values
protein_name = df.iloc[X_test.index[test_idx]]['target_name']

print(f"\n🎯 Generating aptamers for: {protein_name[:60]}...")

# Generate candidates
print(f"\nGenerating 10 candidate DNA aptamers...")
candidates = generate_aptamer_candidates(
    protein_features,
    model,
    kmer_library_3,
    transition_matrix_3,
    output_features,
    n_candidates=10,
    k=3
)

# Display top 5
print(f"\n🏆 TOP 5 GENERATED DNA APTAMERS:")
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
    print(f"   Match Errors: Length±{val['matches']['length_diff']:.1f}bp, "
          f"GC±{val['matches']['gc_diff']:.3f}")

# Compare to real aptamer if available
real_apt_idx = df.index[X_test.index[test_idx]]
if pd.notna(df.loc[real_apt_idx, 'aptamer_sequence']):
    real_apt = df.loc[real_apt_idx, 'aptamer_sequence']
    print(f"\n📋 ACTUAL DNA APTAMER FOR THIS PROTEIN:")
    print(f"   Sequence: {real_apt}")
    print(f"   Length: {len(real_apt)} bp")
    real_gc = (real_apt.count('G') + real_apt.count('C')) / len(real_apt)
    print(f"   GC Content: {real_gc:.2%}")

print("\n" + "="*80)
print("✅ COMPLETE! DNA APTAMER GENERATION SUCCESSFUL")
print("="*80)
print("\n💡 All generated sequences are DNA (A, C, G, T only)")
print("   No RNA/DNA mixing issues!")
print("\n📁 Model saved as: my_aptamer_predictor.pkl")
print("   Use this for generating more aptamers later!")
