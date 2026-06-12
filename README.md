## Sequence Analysis and Phylogenetics

A Python-based bioinformatics project for protein sequence retrieval, FASTA processing, homolog identification, multiple sequence alignment, evolutionary distance calculation, phylogenetic tree construction, sugar-metabolism annotation, and ancestral active-site reconstruction.
This repository is based on a Jupyter Notebook workflow for Sequence Analysis and Phylogenetics.


## Project Overview

This project analyses biological protein sequences using Python. The workflow begins with retrieving protein sequences from UniProt, processing FASTA files, identifying homologous sequences, aligning protein sequences, calculating evolutionary distances, and constructing phylogenetic trees.
The final part focuses on the MALS protein family in yeasts. MALS proteins are involved in carbohydrate metabolism, especially sugars such as maltose and isomaltose. The project investigates how sequence variation, active-site changes, and phylogenetic relationships may explain functional differences among MALS proteins.

## Main Objectives

Retrieve protein sequences from online databases.
Process and interpret FASTA formatted sequence files.
Identify homologous proteins using sequence similarity.
Generate and analyse multiple sequence alignments.
Calculate evolutionary distances using fractional, Poisson, and Gamma correction methods.
Build phylogenetic trees using UPGMA.
Annotate trees with sugar metabolism data.
Infer ancestral active-site residues using maximum parsimony.



---
Software Requirements
This project uses Python and Jupyter Notebook.
Required Python packages:
```text
numpy
matplotlib
jupyter
biopython
```
Install packages inside Jupyter or terminal using:
```python
pip install numpy matplotlib jupyter biopython
```
---
How to Run the Project
Open Jupyter Notebook and run:
```python
from guide import *
import csv
import numpy as np
import matplotlib.pyplot as plt
```
Then open:
```text
workshop_2_assesment.ipynb
```
## Run each section in order.
---
**Project Workflow**

Part A: Database Access

Exercise 1: Retrieve and Process a Single Protein
This exercise retrieves a single protein sequence from UniProt and extracts useful information from the annotation.
Main tasks:
Retrieve a protein sequence.
Print sequence name and annotation.
Identify useful features such as gene name, species, or signal peptide.
Example code:
```python
seq = getSequence('Q9Y2P8', debug=False)

print(seq.name)
print(seq.annot)
print(seq.sequence[:100])
```
---
Exercise 2: Retrieve Multiple Proteins
This exercise searches UniProt for multiple proteins from the same family.
Example code:
```python
seqs = []

names = searchSequences('(reviewed:true) AND (protein_name:DAO)')

for name in names:
    seq = getSequence(name, debug=False)
    seqs.append(seq)
    print(seq.name, seq.annot)

writeFastaFile('dao.fa', seqs)
```
This saves the retrieved DAO protein sequences into a FASTA file.

---
Exercise 3: Processing FASTA Format
FASTA is a simple biological sequence format. Each sequence begins with a header line starting with `>`, followed by one or more sequence lines.
Example:
```fasta
>Protein_1 OS=Homo sapiens
MKRHYEAVVIGGGIIGSAIAYYLAKENK
```
A FASTA file may fail to load if:
The sequence contains invalid characters.
Header lines do not start with `>`.
Sequence formatting is incorrect.
The selected alphabet does not match the sequence type.

---
Part B: Sequence Alignment

Exercise 4: Find Two Putative Homologs
This exercise compares one DAO protein sequence with all other DAO sequences and finds the two closest matches.
Example code:
```python
best1 = 0
idx1 = 0
best2 = 0
idx2 = 0

for idx in range(len(seqs)):
    if idx != my_sequence_number:
        seq = seqs[idx]
        aln = align(seqs[my_sequence_number], seq, blosum62, -7)
        percent = scoreAlignment(aln) / aln.alignlen

        if best1 < percent:
            best2 = best1
            idx2 = idx1
            best1 = percent
            idx1 = idx

        elif best2 < percent:
            best2 = percent
            idx2 = idx

print('Closest homolog:', seqs[idx1].name, best1 * 100)
print('Second closest homolog:', seqs[idx2].name, best2 * 100)
```

Interpretation:
The two sequences with the highest identity values are the closest putative homologs.

---
Exercise 5: Consensus Sequence
A consensus sequence represents the most common amino acid at each alignment column.
```python
def getConsensusForColumn(aln, colidx):
    symcnt = {}
    for seq in aln.seqs:
        mysym = seq[colidx]
        symcnt[mysym] = symcnt.get(mysym, 0) + 1

    consensus = max(symcnt, key=symcnt.get)
    return consensus


def getConsensus(aln):
    consensus_str = ''
    for colidx in range(aln.alignlen):
        consensus_str += getConsensusForColumn(aln, colidx)
    return Sequence(consensus_str, Protein_Alphabet, name='Consensus', gappy=True)
```
---
Exercise 6: Evolutionary Distance Metrics
This exercise adds Poisson and Gamma corrected distances to `calcDistances()` in `guide.py`.
Fractional distance measures the observed proportion of different sites.
Poisson correction accounts for multiple substitutions at the same site:
```python
dist = -math.log(1 - p)
```
Gamma correction accounts for variable rates across sites:
```python
dist = a * ((1 - p) ** (-1.0 / a) - 1)
```
---
Part C: Phylogenetic Analysis

Exercise 7: Select MALS Proteins With Sugar Data
This exercise loads MALS sequences and filters them using species names from `sugars.csv`.
```python
yeasts = dict()

with open('sugars.csv', 'rt') as csvfile:
    reader = csv.reader(csvfile)
    for row in reader:
        yeasts[row[0]] = [y == 'True' for y in row[1:]]

yeast_names = set(yeasts.keys())

seqs = readFastaFile('MalS.fa', Protein_Alphabet)
selected = []

for seq in seqs:
    all_words = [seq.name] + seq.annot.split()
    species = all_words[0]

    if species in yeast_names:
        new_name = '|'.join(all_words)
        new_seq = Sequence(seq.sequence, Protein_Alphabet, name=new_name)
        selected.append(new_seq)

writeFastaFile('select.fa', selected)

print('Selected', len(selected), 'sequences')
```
---
Exercise 8: Consensus of MALS Proteins
After aligning `select.fa` with ClustalX, save the alignment as `select.aln`.
```python
aln_mals = readClustalFile('select.aln', Protein_wGAP)
consensus = getConsensus(aln_mals)
print(consensus)
```
The consensus sequence shows conserved amino acid positions among selected MALS proteins.

---
Exercise 9: Active-Site Mapping
Known active-site columns from the reference MALS alignment are:
```text
173
231-234
294-295
324
437
```
These residues are important for substrate recognition and enzyme activity.
Example answer:
```text
The active-site residues were mapped from the reference MalS alignment to the selected alignment. These positions represent residues involved in sugar binding and may explain differences in substrate specificity among MALS proteins.
```
---
Exercise 10: Build UPGMA Phylogenetic Tree
```python
aln_mals = readClustalFile('select.aln', Protein_wGAP)

tree = runUPGMA(aln_mals, 'poisson')

print(tree)

writeNewickFile('select_upgma.nwk', tree)
```
The output is a Newick tree. It can be copied into a tree viewer or displayed using Python.

---
Exercise 11: Annotate Tree With Sugar Metabolism
Each species is labelled as `True` or `False` depending on whether it can metabolise a sugar.
Example interpretation:
```text
Isomaltose metabolism showed a clear phylogenetic pattern. Saccharomyces species clustered together and were mostly labelled True, suggesting that isomaltose metabolism may have evolved in a shared ancestor of this clade.
```
---
Exercise 12: Ancestral Active-Site Reconstruction
This exercise uses maximum parsimony to infer ancestral active-site residues.
```python
aln_gap = readClustalFile('select.aln', Protein_wGAP)

tree12 = readNewickFile('select_upgma.nwk')

tree12.putAlignment(aln_gap)

tree12.parsimony()
```
To print active-site residues at each node:
```python
def strSites(node, site_columns):
    if node.sequence is not None:
        site_str = ''.join([node.sequence[c] for c in site_columns])
    else:
        site_str = '?' * len(site_columns)

    if node.left is None and node.right is None:
        return '%s_%s' % (node.label, site_str)

    left_str = ''
    right_str = ''

    if node.left is not None:
        left_str = strSites(node.left, site_columns)

    if node.right is not None:
        right_str = strSites(node.right, site_columns)

    return '(%s,%s)%s_%s' % (left_str, right_str, node.label, site_str)

active_site_cols = [168, 226, 227, 228, 229, 289, 290, 319, 432]

full_site_newick = strSites(tree12.root, active_site_cols) + ';'

print(full_site_newick)
```
**Interpretation:**
The labels at the leaf nodes show modern active-site residues. The labels at internal nodes represent ancestral active-site residues inferred using maximum parsimony.

---
**Final Biological Interpretation**

The analysis shows that protein sequence variation is connected to evolutionary and functional differences. DAO protein analysis identifies close homologs based on sequence similarity. MALS protein analysis shows that closely related yeast proteins often share similar active-site residues and sugar-metabolism patterns.
The phylogenetic tree suggests that MALS proteins evolved through gene duplication and divergence. Conserved active-site residues indicate functional constraints, while changes in specific active-site positions likely contributed to differences in substrate specificity, including the ability to metabolise isomaltose.

---

Example GitHub Description
Python-based bioinformatics workflow for protein sequence analysis, homolog identification, multiple sequence alignment, evolutionary distance calculation, phylogenetic tree construction, sugar-metabolism annotation, and ancestral active-site reconstruction.

---
Suggested Topics
```text
bioinformatics
sequence-analysis
phylogenetics
python
jupyter-notebook
protein-sequences
multiple-sequence-alignment
evolutionary-distance
upgma
newick
mals-proteins
ancestral-reconstruction
```
---
Author
Mahe Alam

---
License
MIT License
