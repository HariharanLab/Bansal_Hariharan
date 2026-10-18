#pairwise identity matrix
from Bio import AlignIO #imports model from biopython to read and write seq aln files
import itertools
import csv

# Function to calculate % identity while ignoring gaps
def calculate_identity(seq1, seq2):
    matches, total = 0, 0
    for res1, res2 in zip(seq1, seq2):
        if res1 != "-" and res2 != "-":  # Ignore gap positions
            total += 1
            if res1 == res2:
                matches += 1
    return (matches / total * 100) if total > 0 else 0  # Avoid division by zero

# Function to process alignment file and output CSV
def pairwise_identity_matrix(input_file, file_format, output_csv):
    alignment = AlignIO.read(input_file, file_format)
    species = [record.id for record in alignment]

    with open(output_csv, "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow([""] + species)  # Header row

        for seq1 in alignment:
            row = [seq1.id]
            for seq2 in alignment:
                identity = calculate_identity(seq1.seq, seq2.seq)
                row.append(f"{identity:.2f}")
            writer.writerow(row)

    print(f"Pairwise identity matrix saved to {output_csv}")

# Example usage
    input_file = "alignment.fasta"  # Change to your filename
    file_format = "fasta"  # Change to "clustal" for .aln or "phylip" for .phy
    output_csv = "pairwise_identity.csv"

#input_file = "/Users/twinklebansal/Downloads/LATS2_Gblocks_Cleaned_sequences.fasta"  # Change to your filename
#file_format = "fasta"  # Change to "clustal" for .aln or "phylip" for .phy
#output_csv = "LATS2_clustal_gblocks_aln.csv"
input_file = "/Users/twinklebansal/Downloads/h_c_a_m.aln"  # Change to your filename
file_format = "fasta"  # Change to "clustal" for .aln or "phylip" for .phy
output_csv = "human_agouti_cap_mouse_aln.csv"
pairwise_identity_matrix(input_file, file_format, output_csv)
