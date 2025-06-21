import os
import sys
import subprocess
from Bio import SeqIO, AlignIO
from Bio.Seq import Seq

if len(sys.argv) < 2:
    print("Usage: python3 main_pipeline.py <input_fasta>")
    sys.exit(1)

INPUT_FASTA = sys.argv[1]
OUTPUT_DIR = "output2"
PROTEIN_DIR = os.path.join(OUTPUT_DIR, "translated_proteins")
MSA_INPUT = os.path.join(OUTPUT_DIR, "all_proteins.fasta")
MSA_OUTPUT = os.path.join(OUTPUT_DIR, "msa_alignment.aln")
MUTATION_FILE = os.path.join(OUTPUT_DIR, "mutation_report.txt")
SIFT_FILE = os.path.join(OUTPUT_DIR, "sift_input.txt")

os.makedirs(PROTEIN_DIR, exist_ok=True)
os.makedirs(OUTPUT_DIR, exist_ok=True)

def translate_and_save(record):
    protein = record.seq.translate(to_stop=False)
    out_file = os.path.join(PROTEIN_DIR, f"{record.id}_protein.fasta")
    with open(out_file, "w") as f:
        f.write(f">{record.id}\n{protein}\n")
    return record.id, str(protein)

def write_all_proteins(protein_seqs):
    with open(MSA_INPUT, "w") as f:
        for id, prot in protein_seqs:
            f.write(f">{id}\n{prot}\n")

def run_muscle(infile, outfile):
    print(" Running MUSCLE for MSA...")
    result = subprocess.run(["muscle", "-in", infile, "-out", outfile],
                            capture_output=True, text=True)
    if result.returncode != 0:
        print(" MUSCLE failed:\n", result.stderr)
        sys.exit(1)
    print(f" MSA saved to {outfile}")

def detect_mutations(aln_file):
    alignment = AlignIO.read(aln_file, "fasta")
    ref_seq = alignment[0].seq
    ref_id = alignment[0].id
    mutations = []

    for record in alignment[1:]:
        for i, (ref_aa, alt_aa) in enumerate(zip(ref_seq, record.seq), 1):
            if ref_aa != alt_aa and ref_aa != "-" and alt_aa != "-":
                mutations.append(f"{record.id}: {ref_aa}{i}{alt_aa}")

    with open(MUTATION_FILE, "w") as f:
        f.write("\n".join(mutations))
    with open(SIFT_FILE, "w") as f:
        f.write("\n".join(mutations))
    print(f" Detected mutations saved to {MUTATION_FILE}")
    return mutations

def main():
    protein_seqs = []
    print(" Translating sequences...")
    for record in SeqIO.parse(INPUT_FASTA, "fasta"):
        protein_seqs.append(translate_and_save(record))

    write_all_proteins(protein_seqs)
    run_muscle(MSA_INPUT, MSA_OUTPUT)
    detect_mutations(MSA_OUTPUT)

    print("\n Pipeline complete!")



if __name__ == "__main__":
    main()
