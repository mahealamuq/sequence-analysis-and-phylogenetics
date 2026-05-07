# MALS-Sequence-Analysis-and-Phylogenetics

This repository contains a bioinformatics workflow for analysing the MALS protein family in Ascomycota. MALS genes encode α-glucosidases involved in carbohydrate metabolism, especially sugars such as maltose and isomaltose.  

Complete bioinformatics workflow for MALS protein homolog identification, conserved domain analysis, active-site prediction, multiple sequence alignment, and phylogenetic tree construction using BLASTP, Pfam, MUSCLE, WebLogo, and MEGA in Ascomycota fungal species.

This repository contains a complete bioinformatics workflow for analysing the evolutionary relationships and functional conservation of MALS proteins in Ascomycota fungi. The pipeline identifies homologous protein sequences, detects conserved domains, analyses active-site residues, performs multiple sequence alignment, and constructs bootstrap-supported phylogenetic trees.

The workflow is designed for:
- Protein sequence analysis
- Comparative genomics
- Conserved domain identification
- Multiple sequence alignment
- Active-site analysis
- Molecular evolution studies
- Phylogenetic analysis

## Table of contents

- [Introduction](#introduction)
- [Workflow Overview](#workflow-overview)
- [Software Requirements](#software-requirements)
- [Installation](#installation)
- [Data Download](#data-download)
- [Running the BLASTP Analysis](#running-the-blastp-analysis)
- [Pfam Domain Analysis](#pfam-domain-analysis)
- [Multiple Sequence Alignment](#multiple-sequence-alignment)
- [Sequence Logo Analysis](#sequence-logo-analysis)
- [Active-Site Identification](#active-site-identification)
- [Phylogenetic Tree Construction](#phylogenetic-tree-construction)
- [Results](#results)
- [References](#references)

# MALS Sequence Analysis & Phylogenetics Pipeline

## Introduction

This repository contains a complete sequence analysis and phylogenetics workflow for studying MALS protein evolution and functional divergence in fungal species.

The project focuses on:
- MALS protein family
- Ascomycota fungi
- Conserved catalytic residues
- Protein evolution
- Active-site variation
- Phylogenetic relationships

Protein sequence analysis is used to identify:
- Homologous proteins
- Conserved amino acid regions
- Functional domains
- Catalytic motifs
- Evolutionary conservation

Phylogenetic analysis is used to identify:
- Evolutionary relationships
- Gene duplication events
- Paralogs and orthologs
- Functional divergence
- Bootstrap-supported evolutionary clusters

This workflow demonstrates:
- BLASTP homolog identification
- Pfam conserved domain analysis
- MUSCLE multiple sequence alignment
- WebLogo motif analysis
- Active-site prediction
- Maximum Likelihood phylogenetic analysis

The analyses include:
- Conserved domain identification
- Active-site classification
- Bootstrap phylogenetic trees
- Sequence conservation analysis
- Evolutionary interpretation

---

## Workflow Overview

### Homolog Identification Workflow

```text
Query MALS Protein Sequence
            ↓
BLASTP Search Against UniProtKB/Swiss-Prot
            ↓
Identification of Ascomycota Homologs
            ↓
Selection of Outgroup Sequence
```

---

### Conserved Domain & Alignment Workflow

```text
Protein Sequences
        ↓
Pfam Conserved Domain Search
        ↓
MUSCLE Multiple Sequence Alignment
        ↓
Conserved Residue Analysis
        ↓
WebLogo Sequence Motif Analysis
```

---

### Phylogenetic Analysis Workflow

```text
Aligned Protein Sequences
            ↓
MEGA Phylogenetic Analysis
            ↓
Maximum Likelihood Tree Construction
            ↓
Bootstrap Analysis (100 Replicates)
            ↓
Evolutionary Interpretation
```

---

## Software Requirements

**Operating System**

- Ubuntu 20.04+ recommended

**Required Software**

| Software | Purpose |
|---|---|
| BLASTP | Homolog identification |
| UniProtKB/Swiss-Prot | Curated protein database |
| Pfam | Conserved domain analysis |
| MUSCLE | Multiple sequence alignment |
| MEGA | Phylogenetic analysis |
| WebLogo 3 | Conserved motif visualization |
| SeqKit | FASTA sequence processing |
| SAMtools | Sequence utilities |

---

## Installation

### Update Ubuntu

```bash
sudo apt update
```

---

### Install Required Tools

```bash
sudo apt install -y \
ncbi-blast+ \
muscle \
seqkit \
samtools \
wget \
curl \
unzip
```

---

### Install MEGA

Download from:

```text
https://www.megasoftware.net/
```

---

### Install WebLogo

```bash
pip install weblogo
```

---

## Data Download

### Create Project Directory

```bash
mkdir -p mals_phylogenetics/{data,blast_results,pfam,alignment,weblogo,phylogeny,results,scripts}

cd mals_phylogenetics
```

---

### Download Query Protein Sequence

Example FASTA format:

```fasta
>Query_MALS_Protein
MSTNPKPQRKTKRNTNRRPQDVKFPGGGQIVGGV...
```

Save as:

```text
data/query_sequence.fasta
```

---

### Download Homologous Sequences Using BLASTP

```bash
blastp \
-query data/query_sequence.fasta \
-db swissprot \
-out blast_results/mals_blast_results.txt \
-outfmt 6 \
-remote
```

---

### Download Selected Homolog FASTA Sequences

Save homolog sequences into:

```text
data/ascomycota_homologs.fasta
```

---

## Running the BLASTP Analysis

### Step 1: Run BLASTP

```bash
blastp \
-query data/query_sequence.fasta \
-db swissprot \
-out blast_results/blast_output.txt \
-evalue 1e-5 \
-max_target_seqs 100 \
-outfmt 6 \
-remote
```

---

### Step 2: View BLAST Results

```bash
less blast_results/blast_output.txt
```

---

### Step 3: Select Homologs

Select:
- Ascomycota homologs
- One non-Ascomycota outgroup

Combine into:

```text
data/combined_sequences.fasta
```

---

## Pfam Domain Analysis

### Upload Query Sequence to Pfam

Website:

```text
https://pfam.xfam.org/
```

---

### Record Important Information

Collect:
- Domain name
- Pfam accession
- Clan
- Species distribution
- Available PDB structures

---

### Save Results

```text
pfam/pfam_results.md
```

---

## Multiple Sequence Alignment

### Step 1: MUSCLE Default Alignment

```bash
muscle \
-in data/combined_sequences.fasta \
-out alignment/muscle_default_alignment.fasta
```

---

### Step 2: MUSCLE Stringent Alignment

```bash
muscle \
-gapopen -2.9 \
-gapextend -1.45 \
-in data/combined_sequences.fasta \
-out alignment/muscle_stringent_alignment.fasta
```

---

### Step 3: Compare Alignments

Compare:
- Alignment length
- Conserved sites
- Invariant sites
- Singleton sites
- Parsimony informative sites

---

## Sequence Logo Analysis

### Extract Conserved Domain Region

```bash
seqkit subseq \
-r 1:25 \
alignment/muscle_default_alignment.fasta \
-o weblogo/domain_region.fasta
```

---

### Generate Sequence Logo

```bash
weblogo \
-f weblogo/domain_region.fasta \
-o weblogo/conserved_logo.png
```

---

## Active-Site Identification

### Analyse Conserved Residues

Identify catalytic residues corresponding to:

| Active Site | Predicted Substrate |
|---|---|
| GVA | Maltose |
| SMQ | Isomaltose |

---

### Compare Query and Homologs

Analyse:
- Conserved amino acids
- Catalytic substitutions
- Functional divergence

---

## Phylogenetic Tree Construction

### Step 1: Import Alignment into MEGA

```text
MEGA → Open Alignment File
```

---

### Step 2: Construct Phylogenetic Tree

Select:
- Maximum Likelihood method
- Amino acid substitution model
- Bootstrap method
- 100 bootstrap replications

---

### Step 3: Root Tree Using Outgroup

Use:
- Non-Ascomycota homolog

---

### Step 4: Export Tree

Save:
- Bootstrap consensus tree
- Newick tree file
- Tree image

---

## Results

Expected outputs include:

### BLASTP Results
- Homologous protein list
- Species identification
- Similarity scores
- E-values

### Pfam Results
- Conserved domain identification
- Domain architecture
- Protein family classification

### Alignment Results
- Conserved amino acid regions
- Functional motifs
- Catalytic residues

### Phylogenetic Results
- Maximum Likelihood tree
- Bootstrap support values
- Evolutionary clusters
- Paralog relationships

### Functional Interpretation
- Active-site prediction
- Substrate specificity
- Functional divergence
- Evolutionary conservation

---

## Repository Structure

```text
mals-sequence-analysis-phylogenetics/
│
├── data/
├── blast_results/
├── pfam/
├── alignment/
├── weblogo/
├── phylogeny/
├── results/
├── scripts/
│
├── README.md
└── pipeline.sh
```

---

## Example Output Files

| File | Description |
|---|---|
| blast_output.txt | BLASTP homolog results |
| pfam_results.md | Conserved domain results |
| muscle_default_alignment.fasta | MUSCLE default alignment |
| muscle_stringent_alignment.fasta | Stringent alignment |
| conserved_logo.png | WebLogo motif |
| bootstrap_tree.png | Phylogenetic tree |

---

## Biological Interpretation

The MALS protein family demonstrates:
- Gene duplication
- Functional diversification
- Evolutionary conservation
- Active-site variation

Closely related proteins may share:
- Similar substrate specificity
- Conserved catalytic residues
- Shared evolutionary ancestry

Longer phylogenetic branches may indicate:
- Faster evolutionary rates
- Greater sequence divergence
- Functional adaptation

---

## References

### BLASTP

Altschul SF et al. (1997).  
Gapped BLAST and PSI-BLAST: a new generation of protein database search programs.  
Nucleic Acids Research.

---

### MUSCLE

Edgar RC. (2004).  
MUSCLE: multiple sequence alignment with high accuracy and high throughput.  
Nucleic Acids Research.

---

### Pfam

Mistry J et al. (2021).  
Pfam: The protein families database in 2021.  
Nucleic Acids Research.

---

### MEGA

Tamura K et al. (2021).  
MEGA11: Molecular Evolutionary Genetics Analysis version 11.  
Molecular Biology and Evolution.

---

### WebLogo

Crooks GE et al. (2004).  
WebLogo: a sequence logo generator.  
Genome Research.

---

## Author

Mahe Alam

---

## License

MIT License
