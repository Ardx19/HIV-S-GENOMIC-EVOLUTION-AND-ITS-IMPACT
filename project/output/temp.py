from Bio import AlignIO
import numpy as np
import matplotlib.pyplot as plt
import math
from scipy.ndimage import uniform_filter1d
from scipy.stats import zscore

def plot_entropy_profile(filepath, file_format="fasta", smoothing_window=10, z_threshold=1.5):
    alignment = AlignIO.read(filepath, file_format)

    def shannon_entropy(column):
        unique_residues = set(column)
        total = len(column)
        entropy = 0.0
        for residue in unique_residues:
            p = column.count(residue) / total
            entropy -= p * math.log2(p)
        return entropy

    raw_entropy = [shannon_entropy([record.seq[i] for record in alignment]) 
                   for i in range(alignment.get_alignment_length())]

    # Smooth entropy
    smoothed_entropy = uniform_filter1d(raw_entropy, size=smoothing_window)

    # Z-score normalization
    z_scores = zscore(smoothed_entropy)

    # Plot
    plt.figure(figsize=(16, 6))
    plt.plot(smoothed_entropy, label="Smoothed Entropy", color='darkblue')
    plt.axhline(np.mean(smoothed_entropy), color='gray', linestyle='--', label="Mean")
    
    # Highlight positions above threshold
    high_variability = [i for i, z in enumerate(z_scores) if z > z_threshold]
    plt.scatter(high_variability, [smoothed_entropy[i] for i in high_variability],
                color='red', label=f"Z > {z_threshold}", s=15)

    plt.title("Shannon Entropy Profile (Smoothed + Normalized)")
    plt.xlabel("Amino Acid Position")
    plt.ylabel("Entropy")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.show()
    plt.plot(z_scores, label="Z-scores", color='orange')
    plt.show()
plot_entropy_profile("msa_alignment.aln", file_format="fasta")
plt.savefig("entropy_profile.png")
plt.close()
