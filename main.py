import requests
import time
from Bio import Entrez, SeqIO, Seq
from Bio.Blast import NCBIWWW, NCBIXML

Entrez.email = "meow@iiitd.ac.in"

accession_data = {
    "PDC1": ("NC_001144.5", "pdc1_sequence.fasta", 232389, 234081),
    "PDC5": ("NC_001144.5", "pdc5_sequence.fasta", 410722, 412414),
    "PDC6": ("NC_001139.9", "pdc6_sequence.fasta", 651289, 652981)
}

def downloadseq(accession_info):
    try:
        accession, filename, start, end = accession_info
        handle = Entrez.efetch(db="nuccore", id=accession, rettype="fasta", retmode="text", seq_start=start, seq_stop=end)
        seq_record = SeqIO.read(handle, "fasta")
        SeqIO.write(seq_record, filename, "fasta")
        handle.close()
        print(f"Sliced sequence saved to {filename}")
    except Exception as e:
        print(f"Error occurred during downloading sequence for {accession_info[0]}: {e}")

def blastcom(query_filename, subject_filename):
    try:
        result_handle = NCBIWWW.qblast("blastn", "nt", open(query_filename).read(), alignments=5)
        blast_records = NCBIXML.parse(result_handle)
        totsnp = 0
        for blast_record in blast_records:
            for alignment in blast_record.alignments:
                for hsp in alignment.hsps:
                    totsnp += hsp.align_length - hsp.identities
        return totsnp
    except Exception as e:
        print(f"Error occurred during BLAST comparison: {e}")
        return 0

def translation(dna_sequence):
    try:
        dna_sequence = ''.join(c for c in dna_sequence if c in 'ATGC')
        dna_sequence = dna_sequence[:len(dna_sequence) // 3 * 3]
        amino_acid_sequence = Seq.Seq(dna_sequence).translate()
        return str(amino_acid_sequence)
    except Exception as e:
        print(f"Error occurred during translation: {e}")
        return ""

def polyphen(amino_acid_sequence):
    try:
        url = "https://genetics.bwh.harvard.edu/pph2/"
        data = {
            "seq": amino_acid_sequence,
            "submit": "Run PolyPhen-2"
        }
        response = requests.post(url, data=data, timeout=30)
        if response.ok:
            return response.text
        else:
            return "Error: PolyPhen prediction failed"
    except Exception as e:
        print(f"Error occurred during PolyPhen prediction: {e}")
        return "Error: PolyPhen prediction failed"

for gene, accession_info in accession_data.items():
    downloadseq(accession_info)
    totsnp = blastcom(accession_info[1], accession_info[1])
    print(f"Total SNPs found for {gene}: {totsnp}")
    with open(accession_info[1], "r") as f:
        dna_sequence = f.read().strip()
    amino_acid_sequence = translation(dna_sequence)
    print(f"Amino Acid Sequence for {gene}: {amino_acid_sequence}")
    polyphen_prediction = polyphen(amino_acid_sequence)
    print(f"PolyPhen Prediction for {gene}: {polyphen_prediction}")
