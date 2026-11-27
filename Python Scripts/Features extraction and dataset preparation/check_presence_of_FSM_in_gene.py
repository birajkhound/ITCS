# Article title: Identification of Thermal Class-Specific Genomic and Proteomic Compositional Signatures in Bacteria Using a Weighted Ensemble Feature Ranking Framework, 
# Author Names: Biraj Khound1, Deb K Boruah1, Pranjal Phukan1, Kalyan Sarma1, Prince Das1, Robita Doley2, Alicia Passah1;
# Affiliation: 
# 1 Department of Diagnostic and Interventional Radiology, All India Institute of Medical Sciences, Changsari, Guwahati, 781101, Assam, India.
# 2 Department of Molecular Biology and Biotechnology, Tezpur University, Napaam, Tezpur, 784028, Assam, India.
# E-mail address of the corresponding author: prsbirajaiimsguwahati@gmail.com


import io
import zipfile
import os
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
    # Dictionaries to keep track of counts per file
    # Format: {file_name: count}
    valid_counts = {}
    invalid_counts = {}

    try:
        with zipfile.ZipFile(zip_file_path, "r") as zip_ref:
            file_names = zip_ref.namelist()

            for file_name in file_names:
                if file_name.endswith("cds_from_genomic.fna"):
                    # Initialize counters for this specific file
                    valid_genes = 0
                    invalid_genes = 0
                    current_gene_seq = []

                    with zip_ref.open(file_name) as txt_file:
                        for line_bytes in txt_file:
                            line = line_bytes.decode("utf-8").strip()

                            if line.startswith(">"):
                                # We hit a new header! Process the PREVIOUS gene sequence first
                                if current_gene_seq:
                                    full_seq = "".join(current_gene_seq)
                                    if len(full_seq) % 3 == 0:
                                        valid_genes += 1
                                    else:
                                        invalid_genes += 1
                                    current_gene_seq = (
                                        []
                                    )  # Reset for the new gene
                            else:
                                # Standard sequence line, gather the pieces
                                if line:
                                    current_gene_seq.append(line)

                        # Don't forget to process the very LAST gene in the file
                        if current_gene_seq:
                            full_seq = "".join(current_gene_seq)
                            if len(full_seq) % 3 == 0:
                                valid_genes += 1
                            else:
                                invalid_genes += 1

                    # Store the final counts for this file
                    result = [valid_genes, invalid_genes] 
                    
        return result

    except zipfile.BadZipfile:
        print(f"Error: {zip_file_path} is not a valid zip file.")
        return {}, {}
    except FileNotFoundError:
        print(f"Error: {zip_file_path} not found.")
        return {}, {}
    
def innitiate_excel():
      wt['A1'] = "SN"
      wt['B1'] = "zip file name"
      wt['C1'] = "no of valid genes"
      wt['D1'] = "no of invalid genes"    

def main():
    directory = "D:\\TEMPURA based Dataset\\complete genome sequences"  # Replace with the actual directory path
    xl_file_name_to_save = 'gene details.xlsx'
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
          wt.append([i, organism_zip_file_name, zip_details[0],zip_details[1],])
    wf.save(xl_file_name_to_save)
    print(xl_file_name_to_save+" is saved succesfully!")     
         
if __name__ == '__main__':
    main()          