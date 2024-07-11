I've developed a Python script to analyze the Pyruvate Decarboxylase (PDC) gene from different Saccharomyces cerevisiae strains. This script performs various tasks, including downloading DNA sequences, comparing sequences using BLAST, translating sequences into amino acids, and predicting the impact of single nucleotide polymorphisms (SNPs) on protein function using PolyPhen.

Installation
To use this script, you need Python installed on your system along with the following libraries:

Biopython
Requests
You can install these libraries using pip:
pip install biopython , requests

Set up your email address in the Entrez.email variable to comply with NCBI's guidelines.
Update the accession_data dictionary with the appropriate NCBI accession numbers and sequence ranges for the PDC genes you want to analyze.