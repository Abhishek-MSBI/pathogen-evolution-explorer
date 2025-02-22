import numpy as np
import pandas as pd

def calculate_gc_content(sequence):
    """Calculate GC content of a sequence"""
    sequence = str(sequence).upper()
    gc_count = sequence.count('G') + sequence.count('C')
    total = len(sequence)
    return (gc_count / total * 100) if total > 0 else 0

def calculate_sequence_stats(sequences):
    """Calculate basic statistics for the input sequences"""
    stats = []
    for seq in sequences:
        stats.append({
            'ID': seq.id,
            'Length': len(seq.seq),
            'GC Content (%)': round(calculate_gc_content(seq.seq), 2),
            'Description': seq.description
        })
    return pd.DataFrame(stats)

def compare_sequences(sequences):
    """Compare sequences and calculate similarity metrics"""
    n = len(sequences)
    similarity_matrix = np.zeros((n, n))
    
    for i in range(n):
        for j in range(n):
            matches = sum(a == b for a, b in zip(str(sequences[i].seq), str(sequences[j].seq)))
            similarity = (matches / len(sequences[i].seq)) * 100
            similarity_matrix[i, j] = round(similarity, 2)
    
    return pd.DataFrame(
        similarity_matrix,
        index=[seq.id for seq in sequences],
        columns=[seq.id for seq in sequences]
    )