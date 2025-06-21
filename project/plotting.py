from Bio import AlignIO
import matplotlib.pyplot as plt
import os
import seaborn as sns
import pandas as pd
import numpy as np
from scipy.stats import zscore

OUTPUT_DIR = "output2"
os.makedirs(OUTPUT_DIR, exist_ok=True)

def plot_mutation_profile(aln_file):
    alignment = AlignIO.read(aln_file, "fasta")
    ref_seq = alignment[0].seq
    num_seqs = len(alignment)
    seq_length = len(ref_seq)
    mutation_counts = [0] * seq_length

    for record in alignment[1:]:
        for i in range(seq_length):
            if ref_seq[i] != record.seq[i] and ref_seq[i] != "-" and record.seq[i] != "-":
                mutation_counts[i] += 1

    # Z-normalization
    z_scores = zscore(mutation_counts)
    positions = range(1, seq_length + 1)

    plt.figure(figsize=(14, 5))
    plt.bar(positions, z_scores, color="darkblue")
    plt.axhline(y=1.5, color='red', linestyle='--', label='Z > 1.5 (High mutation)')
    plt.xlabel("Amino Acid Position")
    plt.ylabel("Z-Score of Mutation Count")
    plt.title("Normalized Mutation Frequency (Z-Score)")
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, "mutation_plot_z.png"))
    plt.close()
    print("📊 Z-normalized mutation profile saved as mutation_plot_z.png")

def plot_mutation_heatmap(aln_file):
    alignment = AlignIO.read(aln_file, "fasta")
    ref_seq = alignment[0].seq
    sequence_ids = [record.id for record in alignment]
    seq_length = len(ref_seq)

    heatmap_data = []

    for record in alignment:
        row = []
        for i in range(seq_length):
            ref_aa = ref_seq[i]
            cur_aa = record.seq[i]
            is_mutation = 0 if (ref_aa == cur_aa or cur_aa == "-" or ref_aa == "-") else 1
            row.append(is_mutation)
        heatmap_data.append(row)

    heatmap_array = np.array(heatmap_data)

    # Normalize across columns (positions)
    normed_heatmap = zscore(heatmap_array, axis=0)
    normed_heatmap = np.nan_to_num(normed_heatmap)  # replace NaN (constant columns) with 0

    plt.figure(figsize=(15, len(sequence_ids)))
    sns.heatmap(normed_heatmap, cmap="coolwarm", center=0, cbar_kws={"label": "Z-Score of Mutation"})
    plt.xlabel("Amino Acid Position")
    plt.ylabel("Strain ID")
    plt.title("Z-Normalized Mutation Heatmap")
    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, "mutation_heatmap_z.png"))
    plt.close()
    print("📊 Z-normalized heatmap saved as mutation_heatmap_z.png")


if __name__ == "__main__":
    plot_mutation_profile("output2/msa_alignment.aln")
    plot_mutation_heatmap("output2/msa_alignment.aln")