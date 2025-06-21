# HIV Drug Resistance Analysis Pipeline

This bioinformatics pipeline analyzes the **HIV gag-pol region** across geographically diverse strains to identify mutations, assess structural impacts, and correlate them with drug resistance. The project explores how **genomic evolution** affects treatment strategies by integrating sequence alignment, mutation analysis, entropy visualization, structural comparisons, and external drug resistance databases.

*A more detailed report, including figures and interpretative analysis, is available in the included presentation PDF (`Genomic Evolution of HIV and Its Impact Presentation.pdf`) within this repository.*

---

## 🐝 Project Overview

Human Immunodeficiency Virus (HIV) evolves rapidly due to high mutation and recombination rates, making treatment and vaccine development highly challenging. This pipeline focuses on:

* **Tracking HIV genomic evolution**
* **Identifying region-specific mutation patterns**
* **Analyzing drug resistance hotspots**
* **Inferring structural and functional consequences of mutations**

We particularly study the **Gag-Pol polyprotein** which plays a vital role in virus replication and is subject to drug-induced selective pressures.

---

## 🌍 Geographic Focus

The analysis is conducted using sequences from 5 countries with significant HIV burdens:
**India, Nigeria, South Africa, Brazil, USA**

---

## 🧪 Workflow Summary

1. **Sequence Translation**
   Translates HIV gag-pol nucleotide sequences into proteins.

2. **Multiple Sequence Alignment (MSA)**
   Aligns all protein sequences to a reference for mutation comparison.

3. **Mutation Detection**
   Finds amino acid mutations and logs them in reportable format.

4. **Entropy and Z-score Analysis**
   Computes site-wise variability to highlight hotspots.

5. **Phylogenetic Analysis**
   Constructs evolutionary trees using MEGA to trace common ancestors.

6. **Structural Modeling & RMSD**
   Visualizes protein 3D structures and compares conformations using PyMOL and AlphaFold predictions.

7. **Drug Resistance Prediction**
   Evaluates mutation impact using tools like SIFT and Stanford HIVDB.

---

## 📂 Input and Output

### 🔹 Input

* FASTA file with gag-pol region nucleotide sequences.
  Example: `gag-pol_USA.fasta`

### 🔹 Output (`output2/`)

* `translated_proteins/`: Individual translated sequences.
* `all_proteins.fasta`: Combined protein file for alignment.
* `msa_alignment.aln`: Aligned sequences using MUSCLE.
* `mutation_report.txt`: List of detected mutations.
* `sift_input.txt`: Mutation list formatted for SIFT.
* `mutation_plot_z.png`: Z-score normalized mutation frequency plot.
* `mutation_heatmap_z.png`: Heatmap of mutation Z-scores across strains.

---

## 📂 Project Structure

This repository is organized as follows:

-   **`/` (Root)**
    -   `main_pipeline.py`: The main script that executes the core bioinformatics pipeline, from translation to mutation detection.
    -   `plotting.py`: A script dedicated to generating visualizations from the pipeline's results.
    -   `requirements.txt`: Lists all the Python packages required for the project.
    -   `README.md`: This documentation file.
    -   `hiv_sequences.fasta`: An example FASTA file containing HIV nucleotide sequences.

-   **`/gag-pol/`**
    -   Contains the primary input data, with FASTA files of HIV gag-pol sequences for different geographic regions (e.g., `gag-pol_USA.fasta`).

-   **`/output/` & `/translated_proteins/`**
    -   These directories store the pipeline's outputs. `translated_proteins/` holds the amino acid sequences, while `output/` contains alignments, mutation reports, and plots.

-   **`/sift/`**
    -   Used for storing files and results related to SIFT analysis, which predicts the functional impact of amino acid substitutions.

-   **`/DrugResistanceSdb/`**
    -   Holds data and reports related to the Stanford HIV Drug Resistance Database analysis.

-   **`/pgmol/`**
    -   Contains files for 3D structural analysis, likely including PyMOL (`.pml`, `.pse`) scripts or protein structures (`.pdb`).

---

## 📊 Key Results

* **Entropy peaks** were observed at amino acid regions \~20–40, 80–100, 120–140, and 480–510, indicating potential mutation hotspots.
* **Phylogenetic analysis** revealed evolutionary closeness between India and Brazil, and between Nigeria and the USA.
* **RMSD values** from PyMOL show structural similarities (low RMSD) and differences (high RMSD) with functional implications for drug resistance.
* **Integrase domain conservation** between India and Brazil suggests shared susceptibility to Integrase Strand Transfer Inhibitors (INSTIs).

---

## 🔧 Installation and Usage

### Install Python dependencies

```bash
pip install -r requirements.txt
```

### Ensure MUSCLE is installed

Download and configure MUSCLE from [here](https://www.drive5.com/muscle/manual/install.html).

### Run the main pipeline

```bash
python main_pipeline.py <your_input_file>.fasta
```

### Generate visualizations

```bash
python plotting.py
```

---

## 🔬 Further Analysis

* **SIFT (Mutation Effect)**:
  Submit `sift_input.txt` at [SIFT Web Tool](https://sift.bii.a-star.edu.sg/www/SIFT_seq_submit2.html)

* **Stanford HIV Drug Resistance Database**:
  Upload mutations or full sequences at [Stanford HIVDB](https://hivdb.stanford.edu/)

* **3D Structural Tools**:

  * AlphaFold for structure prediction: [https://alphafold.ebi.ac.uk/](https://alphafold.ebi.ac.uk/)
  * PyMOL for structure visualization and RMSD: [https://pymol.org/](https://pymol.org/)

* **Phylogenetics**:

  * MEGA X software for evolutionary analysis: [https://www.megasoftware.net/](https://www.megasoftware.net/)

---

## 📌 Goal of the Project

To investigate **how geography affects HIV's genetic evolution** and **drug resistance patterns**, with a focus on:

* Identifying conserved vs variable genomic regions
* Mapping structural consequences of mutations
* Improving region-specific treatment strategies

---

## 📚 References & Resources

* [NCBI GenBank](https://www.ncbi.nlm.nih.gov/genbank/)
* [MEGA X](https://www.megasoftware.net/)
* [SIFT Tool](https://sift.bii.a-star.edu.sg/)
* [AlphaFold](https://alphafold.ebi.ac.uk/)
* [Stanford HIVDB](https://hivdb.stanford.edu/)
* [PyMOL](https://pymol.org/)
