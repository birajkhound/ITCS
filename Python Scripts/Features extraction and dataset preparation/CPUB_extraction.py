
# Article title: Identification of Thermal Class-Specific Genomic and Proteomic Compositional Signatures in Bacteria Using a Weighted Ensemble Feature Ranking Framework, 
# Author Names: Biraj Khound1, Deb K Boruah1, Pranjal Phukan1, Kalyan Sarma1, Prince Das1, Robita Doley2, Alicia Passah1;
# Affiliation: 
# 1 Department of Diagnostic and Interventional Radiology, All India Institute of Medical Sciences, Changsari, Guwahati, 781101, Assam, India.
# 2 Department of Molecular Biology and Biotechnology, Tezpur University, Napaam, Tezpur, 784028, Assam, India.
# E-mail address of the corresponding author: prsbirajaiimsguwahati@gmail.com


from collections import defaultdict
import itertools
import os
import zipfile
import numpy as np
import pandas as pd


def get_filenames_in_directory(directory_path):
    try:
        return [f for f in os.listdir(directory_path) if f.endswith(".zip")]
    except FileNotFoundError:
        print(f"Directory not found: {directory_path}")
        return []


def process_gene_sequence(seq_lines, stop_codons):
    """Helper function to clean, check, and repair a single gene sequence."""
    gene_str = "".join(seq_lines).upper().replace("U", "T")

    remainder = len(gene_str) % 3
    if remainder != 0:
        found_stop = False
        # Scan by triplets to find an in-frame stop codon
        for i in range(0, len(gene_str), 3):
            codon = gene_str[i : i + 3]
            if codon in stop_codons:
                new_seq = gene_str[:i]
                if len(new_seq) % 3 == 0:
                    gene_str = new_seq
                    found_stop = True
                    break
        # Fallback: drop the remainder
        if not found_stop:
            gene_str = gene_str[:-remainder]

    return gene_str


def extract_clean_genes(zip_file_path):
    """Extracts CDS sequences and returns a LIST of clean gene strings."""
    try:
        with zipfile.ZipFile(zip_file_path, "r") as zip_ref:
            file_names = zip_ref.namelist()
            for file_name in file_names:
                if file_name.endswith("cds_from_genomic.fna"):
                    valid_master_genes = []
                    current_gene_seq = []
                    stop_codons = {"TAA", "TAG", "TGA"}

                    with zip_ref.open(file_name) as txt_file:
                        for line_bytes in txt_file:
                            line = line_bytes.decode("utf-8").strip()
                            if line.startswith(">"):
                                if current_gene_seq:
                                    clean_gene = process_gene_sequence(
                                        current_gene_seq, stop_codons
                                    )
                                    if clean_gene:
                                        valid_master_genes.append(clean_gene)
                                current_gene_seq = []
                            else:
                                if line:
                                    current_gene_seq.append(line)

                        # Process the final gene in the file
                        if current_gene_seq:
                            clean_gene = process_gene_sequence(
                                current_gene_seq, stop_codons
                            )
                            if clean_gene:
                                valid_master_genes.append(clean_gene)

                    # RETURN THE LIST OF GENES DIRECTLY
                    return valid_master_genes
    except zipfile.BadZipfile:
        print(f"Error: {zip_file_path} is not a valid zip file.")
    except FileNotFoundError:
        print(f"Error: {zip_file_path} not found.")

    return []  # Return empty list if failed


# 1. Standard Genetic Code Setup & Mapping
GENETIC_CODE = {
    "ATA": "I",
    "ATC": "I",
    "ATT": "I",
    "ATG": "M",
    "ACA": "T",
    "ACC": "T",
    "ACG": "T",
    "ACT": "T",
    "AAC": "N",
    "AAT": "N",
    "AAA": "K",
    "AAG": "K",
    "AGC": "S",
    "AGT": "S",
    "AGA": "R",
    "AGG": "R",
    "CTA": "L",
    "CTC": "L",
    "CTG": "L",
    "CTT": "L",
    "CCA": "P",
    "CCC": "P",
    "CCG": "P",
    "CCT": "P",
    "CAC": "H",
    "CAT": "H",
    "CAA": "Q",
    "CAG": "Q",
    "CGA": "R",
    "CGC": "R",
    "CGG": "R",
    "CGT": "R",
    "GTA": "V",
    "GTC": "V",
    "GTG": "V",
    "GTT": "V",
    "GCA": "A",
    "GCC": "A",
    "GCG": "A",
    "GCT": "A",
    "GAC": "D",
    "GAT": "D",
    "GAA": "E",
    "GAG": "E",
    "GGA": "G",
    "GGC": "G",
    "GGG": "G",
    "GGT": "G",
    "TCA": "S",
    "TCC": "S",
    "TCG": "S",
    "TCT": "S",
    "TTC": "F",
    "TTT": "F",
    "TTA": "L",
    "TTG": "L",
    "TAC": "Y",
    "TAT": "Y",
    "TGC": "C",
    "TGT": "C",
    "TGG": "W",
}

STOP_CODONS = {"TAA", "TAG", "TGA"}
bases = ["A", "T", "C", "G"]

# Generate 61 sense codons and 3721 pairs
all_codons = ["".join(c) for c in itertools.product(bases, repeat=3)]
sense_codons = [codon for codon in all_codons if codon not in STOP_CODONS]
valid_pairs = ["".join(p) for p in itertools.product(sense_codons, repeat=2)]

# Map each of the 3721 codon pairs to one of the 400 amino acid pairs
pair_to_aa_group = {}
aa_group_to_pairs = defaultdict(list)

for pair in valid_pairs:
    c1, c2 = pair[:3], pair[3:]
    aa_pair = f"{GENETIC_CODE[c1]}-{GENETIC_CODE[c2]}"
    pair_to_aa_group[pair] = aa_pair
    aa_group_to_pairs[aa_pair].append(pair)

# Set up row headers
headers = ["SN", "zip file name"] + valid_pairs

directory = "D:\\TEMPURA based Dataset\\complete genome sequences"
filenames = get_filenames_in_directory(directory)

result_rows = []

for idx, filename in enumerate(filenames, start=1):
    print(f"Processing {idx} : {filename}")
    file_path = os.path.join(directory, filename)
    genes = extract_clean_genes(file_path)

    if not genes:
        continue

    codon_counts = defaultdict(int)
    pair_counts = defaultdict(int)
    total_pairs_counted = 0

    for gene in genes:
        if len(gene) % 3 != 0 or len(gene) < 6:
            continue

        codons_in_gene = [gene[p : p + 3] for p in range(0, len(gene), 3)]

        for m in range(len(codons_in_gene) - 1):
            c1, c2 = codons_in_gene[m], codons_in_gene[m + 1]
            if (
                c1 in STOP_CODONS
                or c2 in STOP_CODONS
                or c1 not in sense_codons
                or c2 not in sense_codons
            ):
                continue

            codon_counts[c1] += 1
            pair_counts[c1 + c2] += 1
            total_pairs_counted += 1

        if codons_in_gene[-1] in sense_codons:
            codon_counts[codons_in_gene[-1]] += 1

    total_single_codons = sum(codon_counts.values())

    row_data = [idx, filename]

    if total_pairs_counted == 0 or total_single_codons == 0:
        row_data.extend([np.nan] * len(valid_pairs))
        result_rows.append(row_data)
        continue

    # Step A: Pre-calculate baseline expected values for all 3721 pairs
    initial_expected = {}
    for pair in valid_pairs:
        c1, c2 = pair[:3], pair[3:]
        freq_c1 = codon_counts[c1] / total_single_codons
        freq_c2 = codon_counts[c2] / total_single_codons
        initial_expected[pair] = freq_c1 * freq_c2 * total_pairs_counted

    # Step B: Apply the Amino Acid Normalization across the 400 groups
    adjusted_expected = {}
    for aa_group, pairs_in_group in aa_group_to_pairs.items():
        sum_observed = sum(pair_counts[p] for p in pairs_in_group)
        sum_expected = sum(initial_expected[p] for p in pairs_in_group)

        factor = (sum_observed / sum_expected) if sum_expected > 0 else 0.0

        for p in pairs_in_group:
            adjusted_expected[p] = initial_expected[p] * factor

    # Step C: Final Chi-Square pass with low count filter
    for pair in valid_pairs:
        observed = pair_counts[pair]
        expected = adjusted_expected[pair]

        if expected < 3:
            row_data.append(np.nan)
        else:
            chi_sq = ((observed - expected) ** 2) / expected
            row_data.append(round(chi_sq, 6))

    result_rows.append(row_data)

# Export to Excel
df = pd.DataFrame(result_rows, columns=headers)
output_filename = "CPB_values.xlsx"
df.to_excel(output_filename, index=False, engine="openpyxl")
print(f"Success! Balanced file created: {output_filename}")