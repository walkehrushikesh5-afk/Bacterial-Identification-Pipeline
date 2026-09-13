
# 🧬 Bacterial Identification Pipeline Using 16S rRNA Gene Sequences

## 📌 Project Overview

This project presents a bioinformatics pipeline for bacterial identification using 16S rRNA gene sequences.

The pipeline accepts a bacterial DNA sequence in FASTA format and performs sequence analysis followed by similarity-based identification using BLASTN against an NCBI 16S ribosomal RNA reference database.

The application provides both numerical and graphical results through an interactive Streamlit interface.

---

## 🎯 Objective

To develop a simple and automated bioinformatics pipeline that can analyze a bacterial 16S rRNA gene sequence and identify its closest matching bacterial organism based on sequence similarity.

---

## 🔬 Pipeline Workflow

```text
FASTA Sequence
      ↓
FASTA File Reading
      ↓
DNA Sequence Validation
      ↓
Sequence Length
      ↓
A/T/G/C/N Composition
      ↓
GC Percentage
      ↓
Nucleotide Composition Plot
      ↓
BLASTN Search
      ↓
Top 10 Matches
      ↓
Percentage Identity
      ↓
Query Coverage
      ↓
E-value
      ↓
Predicted Organism
      ↓
Identification Assessment
      ↓
Download Results
