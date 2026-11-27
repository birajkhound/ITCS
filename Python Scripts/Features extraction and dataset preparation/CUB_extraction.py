# Article title: Identification of Thermal Class-Specific Genomic and Proteomic Compositional Signatures in Bacteria Using a Weighted Ensemble Feature Ranking Framework, 
# Author Names: Biraj Khound1, Deb K Boruah1, Pranjal Phukan1, Kalyan Sarma1, Prince Das1, Robita Doley2, Alicia Passah1;
# Affiliation: 
# 1 Department of Diagnostic and Interventional Radiology, All India Institute of Medical Sciences, Changsari, Guwahati, 781101, Assam, India.
# 2 Department of Molecular Biology and Biotechnology, Tezpur University, Napaam, Tezpur, 784028, Assam, India.
# E-mail address of the corresponding author: prsbirajaiimsguwahati@gmail.com


import os
import zipfile
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
    genome_details = []
    try:
      with zipfile.ZipFile(zip_file_path, 'r') as zip_ref:
        file_names = zip_ref.namelist()
        # print(f"Files in {zip_file_path}:")
        for file_name in file_names:
          genome_sequence = ""
          no_of_gene = 0
          if file_name.endswith("cds_from_genomic.fna"):      
             with zip_ref.open(file_name) as txt_file:
               for line in txt_file:
                  line = line.decode('utf-8').strip() 
                  if line.startswith('>'):  
                    no_of_gene += 1
                  else:    
                    genome_sequence += line 

               genome_details.extend([no_of_gene,genome_sequence])
               return genome_details  
    except zipfile.BadZipfile:
      print(f"Error: {zip_file_path} is not a valid zip file.")
    except FileNotFoundError:
      print(f"Error: {zip_file_path} not found.")

def access_zip_file(zip_file_path):
    genome_details = []
    try:
        with zipfile.ZipFile(zip_file_path, "r") as zip_ref:
            file_names = zip_ref.namelist()
            for file_name in file_names:
                if file_name.endswith("cds_from_genomic.fna"):
                    no_of_gene = 0
                    valid_master_codons = []
                    current_gene_seq = []
                    stop_codons = {"TAA", "TAG", "TGA"}
                    with zip_ref.open(file_name) as txt_file:
                        for line_bytes in txt_file:
                            line = line_bytes.decode("utf-8").strip()
                            if line.startswith(">"):
                                no_of_gene += 1
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
                    genome_details.extend([no_of_gene,genome_string])
                    return genome_details
    except zipfile.BadZipfile:
        print(f"Error: {zip_file_path} is not a valid zip file.")
    except FileNotFoundError:
        print(f"Error: {zip_file_path} not found.")
    return ""  # Return empty string if no valid file was found or an error occurred
 
 
def calculate_tri_nucleotide_composition(gene_sequence):
    calculated_details = []
    data = gene_sequence.upper()
    data = data.replace('U','T')    
    TTT = 0
    TCT = 0
    TAT = 0
    TGT = 0
    TTC = 0
    TCC = 0
    TAC = 0
    TGC = 0
    TTA = 0
    TCA = 0
    TAA = 0
    TGA = 0
    TTG = 0
    TCG = 0
    TAG = 0
    TGG = 0
    CTT = 0
    CCT = 0
    CAT = 0
    CGT = 0
    CCC = 0
    CAC = 0
    CGC = 0
    CTC = 0
    CCA = 0
    CAA = 0
    CGA = 0
    CTA = 0
    CCG = 0
    CAG = 0
    CGG = 0
    CTG = 0
    ATT = 0
    ACT = 0
    AAT = 0
    AGT = 0
    ACC = 0
    AAC = 0
    AGC = 0
    ATC = 0
    ACA = 0
    AAA = 0
    AGA = 0
    ATA = 0
    ACG = 0
    AAG = 0
    AGG = 0
    ATG = 0
    GTT = 0
    GCT = 0
    GAT = 0
    GGT = 0
    GCC = 0
    GAC = 0
    GGC = 0
    GTC = 0
    GCA = 0
    GAA = 0
    GGA = 0
    GTA = 0
    GCG = 0
    GAG = 0
    GGG = 0
    GTG = 0

    for t in range(0,len(data),3):
         if data[t:t+3] == 'TTT' : 
            TTT += 1 
         elif data[t:t+3] == 'TCT' : 
            TCT += 1
         elif data[t:t+3] == 'TAT' : 
            TAT += 1
         elif data[t:t+3] == 'TGT' : 
            TGT += 1
         elif data[t:t+3] == 'TTC' : 
            TTC += 1
         elif data[t:t+3] == 'TCC' : 
            TCC += 1
         elif data[t:t+3] == 'TAC' : 
            TAC += 1
         elif data[t:t+3] == 'TGC' : 
            TGC += 1
         elif data[t:t+3] == 'TTA' : 
            TTA += 1
         elif data[t:t+3] == 'TCA' : 
            TCA += 1
         elif data[t:t+3] == 'TAA' : 
            TAA += 1
         elif data[t:t+3] == 'TGA' : 
            TGA += 1
         elif data[t:t+3] == 'TTG' : 
            TTG += 1
         elif data[t:t+3] == 'TCG' : 
            TCG += 1
         elif data[t:t+3] == 'TAG' : 
            TAG += 1
         elif data[t:t+3] == 'TGG' : 
            TGG += 1
         elif data[t:t+3] == 'CTT' : 
            CTT += 1
         elif data[t:t+3] == 'CCT' : 
            CCT += 1
         elif data[t:t+3] == 'CAT' : 
            CAT += 1
         elif data[t:t+3] == 'CGT' : 
            CGT += 1
         elif data[t:t+3] == 'CCC' : 
            CCC += 1
         elif data[t:t+3] == 'CAC' : 
            CAC += 1
         elif data[t:t+3] == 'CGC' : 
            CGC += 1
         elif data[t:t+3] == 'CTC' : 
            CTC += 1
         elif data[t:t+3] == 'CCA' : 
            CCA += 1
         elif data[t:t+3] == 'CAA' : 
            CAA += 1
         elif data[t:t+3] == 'CGA' : 
            CGA += 1
         elif data[t:t+3] == 'CTA' : 
            CTA += 1
         elif data[t:t+3] == 'CCG' : 
            CCG += 1
         elif data[t:t+3] == 'CAG' : 
            CAG += 1
         elif data[t:t+3] == 'CGG' : 
            CGG += 1
         elif data[t:t+3] == 'CTG' : 
            CTG += 1 
         elif data[t:t+3] == 'ATT' : 
            ATT += 1
         elif data[t:t+3] == 'ACT' : 
            ACT += 1
         elif data[t:t+3] == 'AAT' : 
            AAT += 1
         elif data[t:t+3] == 'AGT' : 
            AGT += 1
         elif data[t:t+3] == 'ACC' : 
            ACC += 1
         elif data[t:t+3] == 'AAC' : 
            AAC += 1
         elif data[t:t+3] == 'AGC' : 
            AGC += 1
         elif data[t:t+3] == 'ATC' : 
            ATC += 1
         elif data[t:t+3] == 'ACA' : 
            ACA += 1
         elif data[t:t+3] == 'AAA' : 
            AAA += 1
         elif data[t:t+3] == 'AGA' : 
            AGA += 1
         elif data[t:t+3] == 'ATA' : 
            ATA += 1
         elif data[t:t+3] == 'ACG' : 
            ACG += 1
         elif data[t:t+3] == 'AAG' : 
            AAG += 1
         elif data[t:t+3] == 'AGG' : 
            AGG += 1
         elif data[t:t+3] == 'ATG' : 
            ATG += 1
         elif data[t:t+3] == 'GTT' : 
            GTT += 1
         elif data[t:t+3] == 'GCT' : 
            GCT += 1
         elif data[t:t+3] == 'GAT' : 
            GAT += 1
         elif data[t:t+3] == 'GGT' : 
            GGT += 1
         elif data[t:t+3] == 'GCC' : 
            GCC += 1
         elif data[t:t+3] == 'GAC' : 
            GAC += 1
         elif data[t:t+3] == 'GGC' : 
            GGC += 1
         elif data[t:t+3] == 'GTC' : 
            GTC += 1
         elif data[t:t+3] == 'GCA' : 
            GCA += 1
         elif data[t:t+3] == 'GAA' : 
            GAA += 1
         elif data[t:t+3] == 'GGA' : 
            GGA += 1
         elif data[t:t+3] == 'GTA' : 
            GTA += 1 
         elif data[t:t+3] == 'GCG' : 
            GCG += 1
         elif data[t:t+3] == 'GAG' : 
            GAG += 1
         elif data[t:t+3] == 'GGG' : 
            GGG += 1
         elif data[t:t+3] == 'GTG' : 
            GTG += 1
     
    #RSCU_value_calculation
    # RSCU (Relative Synonymous Codon Usage):
    # RSCU = (Observed frequency of a codon) / (Expected frequency of that codon)
    # RSCU values greater than 1 indicate that a codon is used more frequently than expected by chance.   
    # RSCU values less than 1 indicate that a codon is used less frequently than expected.
    if TTT + TTC != 0:
        RSCU_TTT = TTT/((1/2)*(TTT + TTC))
        RSCU_TTC = TTC/((1/2)*(TTT + TTC))
    else:
        RSCU_TTT = 0
        RSCU_TTC = 0    
    # --------------
    if TCT + TCC + TCA + TCG + AGT + AGC != 0:
        RSCU_TCT = TCT/((1/6)*(TCT + TCC + TCA + TCG + AGT + AGC))
        RSCU_TCC = TCC/((1/6)*(TCT + TCC + TCA + TCG + AGT + AGC))
        RSCU_TCA = TCA/((1/6)*(TCT + TCC + TCA + TCG + AGT + AGC))
        RSCU_TCG = TCG/((1/6)*(TCT + TCC + TCA + TCG + AGT + AGC))
        RSCU_AGT = AGT/((1/6)*(TCT + TCC + TCA + TCG + AGT + AGC))
        RSCU_AGC = AGC/((1/6)*(TCT + TCC + TCA + TCG + AGT + AGC))
    else:
        RSCU_TCT = 0
        RSCU_TCC = 0
        RSCU_TCA = 0
        RSCU_TCG = 0
        RSCU_AGT = 0
        RSCU_AGC = 0    
    #--------------
    if TAT + TAC != 0:
        RSCU_TAT = TAT/((1/2)*(TAT + TAC))
        RSCU_TAC = TAC/((1/2)*(TAT + TAC))
    else:
        RSCU_TAT = 0
        RSCU_TAC = 0    
    #-------------
    if TGT + TGC != 0:
        RSCU_TGT = TGT/((1/2)*(TGT + TGC)) 
        RSCU_TGC = TGC/((1/2)*(TGT + TGC))
    else:
        RSCU_TGT = 0
        RSCU_TGC = 0    
    #-------------
    if TTA + TTG + CTT + CTC + CTA + CTG != 0:
        RSCU_TTA = TTA/((1/6)*(TTA + TTG + CTT + CTC + CTA + CTG))
        RSCU_TTG = TTG/((1/6)*(TTA + TTG + CTT + CTC + CTA + CTG))
        RSCU_CTT = CTT/((1/6)*(TTA + TTG + CTT + CTC + CTA + CTG))
        RSCU_CTG = CTG/((1/6)*(TTA + TTG + CTT + CTC + CTA + CTG))
        RSCU_CTC = CTC/((1/6)*(TTA + TTG + CTT + CTC + CTA + CTG))
        RSCU_CTA = CTA/((1/6)*(TTA + TTG + CTT + CTC + CTA + CTG))
    else:
        RSCU_TTA = 0
        RSCU_TTG = 0
        RSCU_CTT = 0
        RSCU_CTG = 0
        RSCU_CTC = 0
        RSCU_CTA = 0    
    #-----------------
    if CCT + CCC + CCA + CCG != 0:
        RSCU_CCT = CCT/((1/4)*(CCT + CCC + CCA + CCG))
        RSCU_CCC = CCC/((1/4)*(CCT + CCC + CCA + CCG))
        RSCU_CCA = CCA/((1/4)*(CCT + CCC + CCA + CCG))
        RSCU_CCG = CCG/((1/4)*(CCT + CCC + CCA + CCG))
    else:
        RSCU_CCT = 0
        RSCU_CCC = 0
        RSCU_CCA = 0
        RSCU_CCG = 0   
    #----------------
    if CAT + CAC != 0:
        RSCU_CAT = CAT/((1/2)*(CAT + CAC))
        RSCU_CAC = CAC/((1/2)*(CAT + CAC))
    else:
        RSCU_CAT = 0
        RSCU_CAC = 0    
    #-------------
    if CGT + CGC + CGA + CGG + AGA + AGG != 0:
        RSCU_CGT = CGT/((1/6)*(CGT + CGC + CGA + CGG + AGA + AGG))
        RSCU_CGC = CGC/((1/6)*(CGT + CGC + CGA + CGG + AGA + AGG))
        RSCU_CGA = CGA/((1/6)*(CGT + CGC + CGA + CGG + AGA + AGG))
        RSCU_CGG = CGG/((1/6)*(CGT + CGC + CGA + CGG + AGA + AGG))
        RSCU_AGA = AGA/((1/6)*(CGT + CGC + CGA + CGG + AGA + AGG))
        RSCU_AGG = AGG/((1/6)*(CGT + CGC + CGA + CGG + AGA + AGG))   
    else:
        RSCU_CGT = 0
        RSCU_CGC = 0
        RSCU_CGA = 0
        RSCU_CGG = 0
        RSCU_AGA = 0
        RSCU_AGG = 0     
    #----------
    if CAA + CAG != 0:
        RSCU_CAA = CAA/((1/2)*(CAA + CAG))
        RSCU_CAG = CAG/((1/2)*(CAA + CAG))
    else:
        RSCU_CAA = 0
        RSCU_CAG = 0    
    #--------
    if AAT + AAC != 0:
        RSCU_AAT = AAT/((1/2)*(AAT + AAC))
        RSCU_AAC = AAC/((1/2)*(AAT + AAC))
    else:
        RSCU_AAT = 0
        RSCU_AAC = 0           
    #-------
    if ATT + ATC + ATA != 0:
        RSCU_ATT = ATT/((1/3)*(ATT + ATC + ATA))
        RSCU_ATC = ATC/((1/3)*(ATT + ATC + ATA))
        RSCU_ATA = ATA/((1/3)*(ATT + ATC + ATA))
    else:
        RSCU_ATT = 0
        RSCU_ATC = 0
        RSCU_ATA = 0    
    #--------
    if ACT + ACC + ACA + ACG != 0:
        RSCU_ACT = ACT/((1/4)*(ACT + ACC + ACA + ACG))
        RSCU_ACA = ACA/((1/4)*(ACT + ACC + ACA + ACG))
        RSCU_ACC = ACC/((1/4)*(ACT + ACC + ACA + ACG))
        RSCU_ACG = ACG/((1/4)*(ACT + ACC + ACA + ACG))
    else:
        RSCU_ACT = 0
        RSCU_ACA = 0
        RSCU_ACC = 0
        RSCU_ACG = 0 
    #-------
    if AAA + AAG != 0:
        RSCU_AAA = AAA/((1/2)*(AAA + AAG))
        RSCU_AAG = AAG/((1/2)*(AAA + AAG))
    else:
        RSCU_AAA = 0
        RSCU_AAG = 0    
    #------
    if GTT + GTC + GTA + GTG != 0:  
        RSCU_GTT = GTT/((1/4)*(GTT + GTC + GTA + GTG))
        RSCU_GTA = GTA/((1/4)*(GTT + GTC + GTA + GTG))
        RSCU_GTC = GTC/((1/4)*(GTT + GTC + GTA + GTG))
        RSCU_GTG = GTG/((1/4)*(GTT + GTC + GTA + GTG))
    else:
        RSCU_GTT = 0
        RSCU_GTA = 0
        RSCU_GTC = 0
        RSCU_GTG = 0    
    #------
    if GTC + GCC + GCA + GCG != 0:
        RSCU_GCT = GCT/((1/4)*(GTC + GCC + GCA + GCG))
        RSCU_GCC = GCC/((1/4)*(GTC + GCC + GCA + GCG))
        RSCU_GCA = GCA/((1/4)*(GTC + GCC + GCA + GCG))
        RSCU_GCG = GCG/((1/4)*(GTC + GCC + GCA + GCG))
    else:
        RSCU_GCT = 0
        RSCU_GCC = 0
        RSCU_GCA = 0
        RSCU_GCG = 0
    #------
    if GAT + GAC != 0:
        RSCU_GAC = GAC/((1/2)*(GAT + GAC))
        RSCU_GAT = GAT/((1/2)*(GAT + GAC))
    else:
        RSCU_GAC = 0
        RSCU_GAT = 0
    #------
    if GGT + GGC + GGA + GGG != 0:
        RSCU_GGT = GGT/((1/4)*(GGT + GGC + GGA + GGG))
        RSCU_GGC = GGC/((1/4)*(GGT + GGC + GGA + GGG))
        RSCU_GGA = GGA/((1/4)*(GGT + GGC + GGA + GGG))    
        RSCU_GGG = GGG/((1/4)*(GGT + GGC + GGA + GGG))
    else:
        RSCU_GGT = 0
        RSCU_GGC = 0
        RSCU_GGA = 0
        RSCU_GGG = 0    
    #-----
    if GAA + GAG != 0:
        RSCU_GAA = GAA/((1/2)*(GAA + GAG))
        RSCU_GAG = GAG/((1/2)*(GAA + GAG))
    else:
        RSCU_GAA = 0
        RSCU_GAG = 0

    calculated_details.extend([ RSCU_TTT ,RSCU_TCT ,RSCU_TAT ,RSCU_TGT ,RSCU_TTC ,RSCU_TCC ,RSCU_TAC ,RSCU_TGC ,RSCU_TTA ,RSCU_TCA ,RSCU_TTG ,RSCU_TCG,RSCU_CTT ,RSCU_CCT ,RSCU_CAT ,RSCU_CGT ,RSCU_CCC ,RSCU_CAC ,RSCU_CGC ,RSCU_CTC ,RSCU_CCA ,RSCU_CAA ,RSCU_CGA ,RSCU_CTA ,RSCU_CCG ,RSCU_CAG ,RSCU_CGG ,RSCU_CTG ,RSCU_ATT ,RSCU_ACT ,RSCU_AAT ,RSCU_AGT ,RSCU_ACC ,RSCU_AAC ,RSCU_AGC ,RSCU_ATC ,RSCU_ACA ,RSCU_AAA ,RSCU_AGA ,RSCU_ATA ,RSCU_ACG ,RSCU_AAG ,RSCU_AGG ,RSCU_GTT ,RSCU_GCT ,RSCU_GAT ,RSCU_GGT ,RSCU_GCC ,RSCU_GAC ,RSCU_GGC ,RSCU_GTC ,RSCU_GCA ,RSCU_GAA ,RSCU_GGA ,RSCU_GTA ,RSCU_GCG ,RSCU_GAG ,RSCU_GGG ,RSCU_GTG ])
    return calculated_details

def innitiate_excel():
      wt['A1'] = "SN"
      wt['B1'] = "zip file name"
      wt['C1'] = "no of genes"
      wt['D1'] = "TTT"
      wt['E1'] = "TCT"
      wt['F1'] = "TAT" 
      wt['G1'] = "TGT"
      wt['H1'] = "TTC"
      wt['I1'] = "TCC"
      wt['J1'] = "TAC"
      wt['K1'] = "TGC"
      wt['L1'] = "TTA"
      wt['M1'] = "TCA"
      wt['N1'] = "TTG"
      wt['O1'] = "TCG"
      wt['P1'] = "CTT"
      wt['Q1'] = "CCT"
      wt['R1'] = "CAT"
      wt['S1'] = "CGT"
      wt['T1'] = "CCC"
      wt['U1'] = "CAC"
      wt['V1'] = "CGC"
      wt['W1'] = "CTC"
      wt['X1'] = "CCA"
      wt['Y1'] = "CAA"
      wt['Z1'] = "CGA"
      wt['AA1'] = "CTA"
      wt['AB1'] = "CCG"
      wt['AC1'] = "CAG"
      wt['AD1'] = "CGG"
      wt['AE1'] = "CTG"
      wt['AF1'] = "ATT"
      wt['AG1'] = "ACT"
      wt['AH1'] = "AAT"
      wt['AI1'] = "AGT"
      wt['AJ1'] = "ACC"
      wt['AK1'] = "AAC"
      wt['AL1'] = "AGC"
      wt['AM1'] = "ATC"
      wt['AN1'] = "ACA"
      wt['AO1'] = "AAA"
      wt['AP1'] = "AGA"
      wt['AQ1'] = "ATA"
      wt['AR1'] = "ACG"
      wt['AS1'] = "AAG"
      wt['AT1'] = "AGG"
      wt['AU1'] = "GTT"
      wt['AV1'] = "GCT"
      wt['AW1'] = "GAT"
      wt['AX1'] = "GGT"
      wt['AY1'] = "GCC"
      wt['AZ1'] = "GAC"
      wt['BA1'] = "GGC"
      wt['BB1'] = "GTC"
      wt['BC1'] = "GCA"
      wt['BD1'] = "GAA"
      wt['BE1'] = "GGA"
      wt['BF1'] = "GTA"
      wt['BG1'] = "GCG"
      wt['BH1'] = "GAG"
      wt['BI1'] = "GGG"
      wt['BJ1'] = "GTG"


def main():
    directory = "D:\\TEMPURA based Dataset\\complete genome sequences"  # Replace with the actual directory path
    xl_file_name_to_save = 'CUB_values.xlsx'
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
            t3 = calculate_tri_nucleotide_composition(zip_details[1]) 
            no_of_genes = zip_details[0]
          else:
            t3 = [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ]
          wt.append([i, organism_zip_file_name, no_of_genes, t3[0], t3[1], t3[2], t3[3], t3[4], t3[5], t3[6], t3[7], t3[8], t3[9], t3[10], t3[11], t3[12], t3[13], t3[14], t3[15], t3[16], t3[17], t3[18], t3[19], t3[20], t3[21], t3[22], t3[23], t3[24], t3[25], t3[26], t3[27], t3[28], t3[29], t3[30], t3[31], t3[32], t3[33], t3[34], t3[35], t3[36], t3[37], t3[38], t3[39], t3[40], t3[41], t3[42], t3[43], t3[44], t3[45], t3[46], t3[47], t3[48], t3[49], t3[50], t3[51], t3[52], t3[53], t3[54], t3[55], t3[56], t3[57], t3[58]])
    wf.save(xl_file_name_to_save)
    print(xl_file_name_to_save+" is saved succesfully!")

if __name__ == '__main__':
    main()

