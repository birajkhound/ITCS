# Article title: Identification of Thermal Class-Specific Genomic and Proteomic Compositional Signatures in Bacteria Using a Weighted Ensemble Feature Ranking Framework, 
# Author Names: Biraj Khound1, Deb K Boruah1, Pranjal Phukan1, Kalyan Sarma1, Prince Das1, Robita Doley2, Alicia Passah1;
# Affiliation: 
# 1 Department of Diagnostic and Interventional Radiology, All India Institute of Medical Sciences, Changsari, Guwahati, 781101, Assam, India.
# 2 Department of Molecular Biology and Biotechnology, Tezpur University, Napaam, Tezpur, 784028, Assam, India.
# E-mail address of the corresponding author: prsbirajaiimsguwahati@gmail.com


import os
import io
from Bio import SeqIO
import dinuq
import zipfile
import tempfile
from openpyxl import Workbook

wf = Workbook()
wt = wf.active

def get_filenames_in_directory(directory_path):
  try:
    return os.listdir(directory_path)
  except FileNotFoundError:
    print(f"Directory not found: {directory_path}")
    return []

def access_zip_file(zip_file_path):
    try:
        with zipfile.ZipFile(zip_file_path, "r") as zip_ref:
            file_names = zip_ref.namelist()
            for file_name in file_names:
                if file_name.endswith("cds_from_genomic.fna"):
                    valid_master_codons = []
                    current_gene_seq = []
                    stop_codons = {"TAA", "TAG", "TGA"}
                    with zip_ref.open(file_name) as txt_file:
                        for line_bytes in txt_file:
                            line = line_bytes.decode("utf-8").strip()
                            if line.startswith(">"):
                                # Process the PREVIOUS gene sequence when a new header is hit
                                if current_gene_seq:
                                    gene_str = "".join(current_gene_seq).upper()
                                    gene_str = gene_str.replace("U", "T")
                                    # Check and repair length if it's not a multiple of 3
                                    remainder = len(gene_str) % 3
                                    if remainder != 0:
                                        gene_str = gene_str[:-remainder]
                                    valid_master_codons.append(gene_str)
                                current_gene_seq = []  # Reset for the new gene
                            else:
                                if line:
                                    current_gene_seq.append(line)
                        # Don't forget to process the very LAST gene in the file
                        if current_gene_seq:
                            gene_str = "".join(current_gene_seq).upper()
                            gene_str = gene_str.replace("U", "T")
                            # --- START OF MODIFIED LOGIC ---
                            remainder = len(gene_str) % 3
                            if remainder != 0:
                                found_stop = False
                                # Scan by triplets to find an in-frame stop codon
                                for i in range(0, len(gene_str), 3):
                                    codon = gene_str[i : i + 3]
                                    if codon in stop_codons:
                                        # Slice right up to (but excluding) this stop codon
                                        new_seq = gene_str[:i]
                                        if len(new_seq) % 3 == 0:
                                            gene_str = new_seq
                                            found_stop = True
                                            break
                                # Fallback: If no stop codon was found or math failed, drop the tail remainder
                                if not found_stop:
                                    gene_str = gene_str[:-remainder]
                            # --- END OF MODIFIED LOGIC ---
                            valid_master_codons.append(gene_str)
                    # Concatenate all corrected genes together into a single string
                    genome_string = "".join(valid_master_codons)
                    # Return immediately since we found and processed the target file
                    return genome_string
    except zipfile.BadZipfile:
        print(f"Error: {zip_file_path} is not a valid zip file.")
    except FileNotFoundError:
        print(f"Error: {zip_file_path} not found.")
    return ""  # Return empty string if no valid file was found or an error occurred

def calculate_di_nucleotide_composition(genome):
    stop_codons = {"TAA", "TAG", "TGA"}
    cleaned_sequence = "".join(
      [
        genome[i : i + 3]
        for i in range(0, len(genome), 3)
        if genome[i : i + 3] not in stop_codons
      ])
    if len(cleaned_sequence)%3 == 0:
        calculated_details = []
        g_count = cleaned_sequence.count('G')  
        c_count = cleaned_sequence.count('C') 
        a_count = cleaned_sequence.count('A') 
        t_count = cleaned_sequence.count('T') 
        length_of_the_seq = len(cleaned_sequence)
        at_skew = (a_count - t_count)/(a_count + t_count)
        gc_skew = (g_count - c_count)/(g_count + c_count)
        gc_percentage = ((g_count + c_count)/(g_count + c_count + a_count + t_count))*100
        at_percentage = ((a_count + t_count)/(g_count + c_count + a_count + t_count))*100  
        calculated_details.extend([length_of_the_seq, a_count, t_count, g_count, c_count, at_skew, gc_skew, at_percentage, gc_percentage])
#-------------------------------------  
        fasta_header = ">my_custom_sequence\n"
        fasta_file = io.StringIO(fasta_header + cleaned_sequence + "TAA")
        alldns = [ 'ApA', 'ApC', 'ApG', 'ApT', 'CpA', 'CpC', 'CpG', 'CpT', 'GpA', 'GpC', 'GpG', 'GpT', 'TpA', 'TpC', 'TpG', 'TpT']
        with tempfile.NamedTemporaryFile(
            mode="w+", suffix=".fasta", delete=False
        ) as temp_fasta:
            temp_fasta.write(fasta_file.getvalue())
            temp_fasta.flush()
            rsduc = dinuq.RSDUc(fasta_file=temp_fasta.name, dinucl=alldns)
            position = ['pos1', 'pos2', 'bridge']
            seq_data = rsduc.get("my_custom_sequence", {})
            for dns in alldns:
                for pos in position:
                    feature_key = dns + pos
                    val_list = seq_data.get(feature_key, [1])  # Fallback to [0.0] if missing
                    calculated_details.extend([val_list[0]])
    return calculated_details
            
def innitiate_excel():
        wt['A1'] = "SN"
        wt['B1'] = "zip file name"
        wt['C1'] = "length of the cds"
        wt['D1'] = "A_(CDS)"
        wt['E1'] = "T_(CDS)"
        wt['F1'] = "G_(CDS)"
        wt['G1'] = "C_(CDS)" 
        wt['H1'] = "AT_SKEW"
        wt['I1'] = "GC_SKEW"
        wt['J1'] = "AT_PERCENTAGE"
        wt['K1'] = "GC_PERCENTAGE"
        wt['L1'] = "ApApos1"
        wt['M1'] = "ApApos2"
        wt['N1'] = "ApAbridge"
        wt['O1'] = "ApCpos1"
        wt['P1'] = "ApCpos2"
        wt['Q1'] = "ApCbridge"
        wt['R1'] = "ApGpos1"
        wt['S1'] = "ApGpos2"
        wt['T1'] = "ApGbridge"
        wt['U1'] = "ApTpos1"
        wt['V1'] = "ApTpos2"
        wt['W1'] = "ApTbridge"
        wt['X1'] = "CpApos1"
        wt['Y1'] = "CpApos2"
        wt['Z1'] = "CpAbridge"
        wt['AA1'] = "CpCpos1"
        wt['AB1'] = "CpCpos2"
        wt['AC1'] = "CpCbridge"
        wt['AD1'] = "CpGpos1"
        wt['AE1'] = "CpGpos2"
        wt['AF1'] = "CpGbridge"
        wt['AG1'] = "CpTpos1"
        wt['AH1'] = "CpTpos2"
        wt['AI1'] = "CpTbridge"
        wt['AJ1'] = "GpApos1"
        wt['AK1'] = "GpApos2"
        wt['AL1'] = "GpAbridge"
        wt['AM1'] = "GpCpos1"
        wt['AN1'] = "GpCpos2"
        wt['AO1'] = "GpCbridge"
        wt['AP1'] = "GpGpos1"
        wt['AQ1'] = "GpGpos2"
        wt['AR1'] = "GpGbridge"
        wt['AS1'] = "GpTpos1"
        wt['AT1'] = "GpTpos2"
        wt['AU1'] = "GpTbridge"
        wt['AV1'] = "TpApos1"
        wt['AW1'] = "TpApos2"
        wt['AX1'] = "TpAbridge"
        wt['AY1'] = "TpCpos1"
        wt['AZ1'] = "TpCpos2"
        wt['BA1'] = "TpCbridge"
        wt['BB1'] = "TpGpos1"
        wt['BC1'] = "TpGpos2"
        wt['BD1'] = "TpGbridge"
        wt['BE1'] = "TpTpos1"
        wt['BF1'] = "TpTpos2"
        wt['BG1'] = "TpTbridge"

def main():
    directory = "D:\\TEMPURA based Dataset\\complete genome sequences"  # Replace with the actual directory path
    xl_file_name_to_save = 'DUB_values.xlsx'
    innitiate_excel()
    filenames = get_filenames_in_directory(directory)
    i = 0
    if filenames:
      for filename in filenames:
          i += 1
          organism_zip_file_name = filename
          print(i, " : ", organism_zip_file_name)
          file_path = os.path.join(directory, filename)
          zip_details = access_zip_file(file_path)
          if zip_details != None:
            t3 = calculate_di_nucleotide_composition(zip_details) 
          else:
            t3 = [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ]
          wt.append([i, organism_zip_file_name, t3[0], t3[1], t3[2], t3[3], t3[4], t3[5], t3[6], t3[7], t3[8], t3[9], t3[10], t3[11], t3[12], t3[13], t3[14], t3[15], t3[16], t3[17], t3[18], t3[19], t3[20], t3[21], t3[22], t3[23], t3[24], t3[25], t3[26], t3[27], t3[28], t3[29], t3[30], t3[31], t3[32], t3[33], t3[34], t3[35], t3[36], t3[37], t3[38], t3[39], t3[40], t3[41], t3[42], t3[43], t3[44], t3[45], t3[46], t3[47], t3[48], t3[49], t3[50], t3[51], t3[52], t3[53], t3[54], t3[55], t3[56]])
    wf.save(xl_file_name_to_save)
    print(xl_file_name_to_save+" is saved succesfully!")

if __name__ == '__main__':
    main()

