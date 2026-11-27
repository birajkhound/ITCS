# Article title: Identification of Thermal Class-Specific Genomic and Proteomic Compositional Signatures in Bacteria Using a Weighted Ensemble Feature Ranking Framework, 
# Author Names: Biraj Khound1, Deb K Boruah1, Pranjal Phukan1, Kalyan Sarma1, Prince Das1, Robita Doley2, Alicia Passah1;
# Affiliation: 
# 1 Department of Diagnostic and Interventional Radiology, All India Institute of Medical Sciences, Changsari, Guwahati, 781101, Assam, India.
# 2 Department of Molecular Biology and Biotechnology, Tezpur University, Napaam, Tezpur, 784028, Assam, India.
# E-mail address of the corresponding author: prsbirajaiimsguwahati@gmail.com


import itertools
import pandas as pd
import os
import zipfile
import math
from collections import Counter

def get_filenames_in_directory(directory_path):
    try:
        return os.listdir(directory_path)
    except FileNotFoundError:
        print(f"Directory not found: {directory_path}")
        return []

def access_zip_file_sequences(zip_file_path):
    """
    Extracts all individual protein sequences from .faa files inside the zip.
    Returns a list of clean protein strings to avoid junction artifacts.
    """
    sequences = []
    try:
        with zipfile.ZipFile(zip_file_path, 'r') as zip_ref:
            for file_name in zip_ref.namelist():
                if file_name.endswith(".faa"):      
                    current_seq = []
                    with zip_ref.open(file_name) as txt_file:
                        for line in txt_file:
                            line = line.decode('utf-8').strip() 
                            if line.startswith('>'):
                                if current_seq:
                                    sequences.append("".join(current_seq).upper())
                                    current_seq = []
                            else:
                                current_seq.append(line)
                        if current_seq:
                            sequences.append("".join(current_seq).upper())
    except zipfile.BadZipfile:
        print(f"Error: {zip_file_path} is not a valid zip file.")
    except FileNotFoundError:
        print(f"Error: {zip_file_path} not found.")
    return sequences

# Define standard amino acids and generate all 400 dipeptides
amino_acids = ['A','G','I','L','P','V','F','W','Y','D','E','R','H','K','S','T','C','M','N','Q']
dipeptide_list = ["".join(p) for p in itertools.product(amino_acids, repeat=2)]

# Setup headers for the final DataFrame
headers = ["SN", "zip file name"] + dipeptide_list

result_rows = []
directory = "D:\\TEMPURA based Dataset\\complete genome sequences"  # Replace with your actual directory path
filenames = get_filenames_in_directory(directory)

if filenames:
    for idx, filename in enumerate(filenames, start=1):
        print(f"{idx} : {filename}")
        file_path = os.path.join(directory, filename)
        
        # Get list of individual protein sequences from this archive
        protein_sequences = access_zip_file_sequences(file_path)
        if not protein_sequences:
            continue

        # Step 1: Count observed occurrences efficiently in a single pass
        observed_counts = Counter()
        first_pos_counts = Counter()
        second_pos_counts = Counter()
        total_dipeptides = 0

        for seq in protein_sequences:
            # Filter sequence to ensure we only look at standard amino acids
            clean_seq = [aa for aa in seq if aa in amino_acids]
            
            # Extract positions and dipeptides
            for j in range(len(clean_seq) - 1):
                pos1 = clean_seq[j]
                pos2 = clean_seq[j+1]
                dipeptide = pos1 + pos2
                
                observed_counts[dipeptide] += 1
                first_pos_counts[pos1] += 1
                second_pos_counts[pos2] += 1
                total_dipeptides += 1

        # Skip files with no valid dipeptides
        if total_dipeptides == 0:
            continue

        # Step 2: Calculate Pi and Pj frequencies
        p_i = {aa: first_pos_counts[aa] / total_dipeptides for aa in amino_acids}
        p_j = {aa: second_pos_counts[aa] / total_dipeptides for aa in amino_acids}

        # Step 3: Compute S_ij for each of the 400 combinations
        row_data = [idx, filename]
        
        for dipeptide in dipeptide_list:
            aa_i = dipeptide[0]
            aa_j = dipeptide[1]            
            # Expected value: E_ij = P_i * P_j * N
            expected = p_i[aa_i] * p_j[aa_j] * total_dipeptides
            observed = observed_counts[dipeptide]            
            # Calculate S_ij normalization (using Vonderviszt standard deviation formula)
            if expected > 0:
                s_ij = observed / expected
            else:
                s_ij = 0.0 # Standard fallback if expected is mathematically zero                
            row_data.append(round(s_ij, 4)) # Rounded for clean presentation in Excel            
        result_rows.append(row_data)

# Create DataFrame and export
df = pd.DataFrame(result_rows, columns=headers)
df.to_excel('AAPB_values.xlsx', index=False, engine='openpyxl')
print("Excel file created with name AAPB_values.xlsx containing Sij values!")