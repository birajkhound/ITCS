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
    protein_sequence = ""
    no_of_proteins = 0
    p_details = []
    try:
      with zipfile.ZipFile(zip_file_path, 'r') as zip_ref:
        file_names = zip_ref.namelist()
        # print(f"Files in {zip_file_path}:")
        for file_name in file_names:
          if file_name.endswith(".faa"):       
             with zip_ref.open(file_name) as txt_file:
               for line in txt_file:
                  line = line.decode('utf-8').strip() 
                  if line.startswith('>'):  
                    no_of_proteins += 1 
                  else:    
                    protein_sequence += line 
             p_details.extend([no_of_proteins,protein_sequence])
      return p_details
    except zipfile.BadZipfile:
      print(f"Error: {zip_file_path} is not a valid zip file.")
    except FileNotFoundError:
      print(f"Error: {zip_file_path} not found.")

def calculator(P_seq):
    p_details = []
    A_no = 0
    G_no = 0
    I_no = 0
    L_no = 0
    P_no = 0
    V_no = 0
    F_no = 0
    W_no = 0
    Y_no = 0
    D_no = 0
    E_no = 0
    R_no = 0
    H_no = 0
    K_no = 0
    S_no = 0
    T_no = 0
    C_no = 0
    M_no = 0
    N_no = 0
    Q_no = 0
    p_seq = P_seq.upper()
    seq_len = len(p_seq)
    for x in range(0,seq_len):
        if p_seq[x:x+1] == 'A' : 
          A_no += 1 
        elif p_seq[x:x+1] == 'G' : 
          G_no += 1
        elif p_seq[x:x+1] == 'I' : 
          I_no += 1
        elif p_seq[x:x+1] == 'L' : 
          L_no += 1
        elif p_seq[x:x+1] == 'P' : 
          P_no += 1
        elif p_seq[x:x+1] == 'V' : 
          V_no += 1
        elif p_seq[x:x+1] == 'F' : 
          F_no += 1
        elif p_seq[x:x+1] == 'W' : 
          W_no += 1
        elif p_seq[x:x+1] == 'Y' : 
          Y_no += 1
        elif p_seq[x:x+1] == 'D' : 
          D_no += 1
        elif p_seq[x:x+1] == 'E' : 
          E_no += 1
        elif p_seq[x:x+1] == 'R' : 
          R_no += 1
        elif p_seq[x:x+1] == 'H' : 
          H_no += 1
        elif p_seq[x:x+1] == 'K' : 
          K_no += 1
        elif p_seq[x:x+1] == 'S' : 
          S_no += 1
        elif p_seq[x:x+1] == 'T' : 
          T_no += 1
        elif p_seq[x:x+1] == 'C' : 
          C_no += 1
        elif p_seq[x:x+1] == 'M' : 
          M_no += 1
        elif p_seq[x:x+1] == 'N' : 
          N_no += 1
        elif p_seq[x:x+1] == 'Q' : 
          Q_no += 1
    known_aa = A_no + G_no + I_no + L_no + P_no + V_no + F_no + W_no + Y_no + D_no + E_no + R_no + H_no + K_no + S_no + T_no + C_no + M_no + N_no + Q_no
    average_aa = known_aa/20
    A = A_no/average_aa
    G = G_no/average_aa
    I = I_no/average_aa
    L = L_no/average_aa
    P = P_no/average_aa
    V = V_no/average_aa
    F = F_no/average_aa
    W = W_no/average_aa
    Y = Y_no/average_aa
    D = D_no/average_aa 
    E = E_no/average_aa
    R = R_no/average_aa
    H = K_no/average_aa
    K = K_no/average_aa
    S = S_no/average_aa
    T = T_no/average_aa
    C = C_no/average_aa
    M = M_no/average_aa
    N = N_no/average_aa
    Q = Q_no/average_aa
    p_details.extend([seq_len, known_aa, A, G, I, L, P, V, F, W, Y, D, E, R, H, K, S, T, C, M, N, Q])
    return p_details

def innitiate_excel():
      wt['A1'] = "SN"
      wt['B1'] = "zip file name"
      wt['C1'] = "no of proteins"
      wt['D1'] = "aa sequence length"
      wt['E1'] = "known aa in the seq"
      wt['F1'] = "A"
      wt['G1'] = "G"
      wt['H1'] = "I"
      wt['I1'] = "L"
      wt['J1'] = "P"
      wt['K1'] = "V"
      wt['L1'] = "F"
      wt['M1'] = "W"
      wt['N1'] = "Y"
      wt['O1'] = "D"
      wt['P1'] = "E"
      wt['Q1'] = "R"
      wt['R1'] = "H"
      wt['S1'] = "K"
      wt['T1'] = "S"
      wt['U1'] = "T"
      wt['V1'] = "C"
      wt['W1'] = "M"
      wt['X1'] = "N"
      wt['Y1'] = "Q"

def main():
    directory = "D:\\TEMPURA based Dataset\\complete genome sequences"  # Replace with the actual directory path
    xl_file_name_to_save = 'AAB_values.xlsx'
    innitiate_excel()
    filenames = get_filenames_in_directory(directory)
    i = 0
    if filenames:
      for filename in filenames:
          i += 1
          organism_zip_file_name = filename
          file_path = os.path.join(directory, filename)
          zip_details = access_zip_file(file_path)
          if len(zip_details) == 2:
             no_of_proteins = zip_details[0]
             t = calculator(zip_details[1]) 
          else:   
            no_of_proteins = 0
            t[0] = 0
            t[1] = 0
            t[2] = 0
            t[3] = 0
            t[4] = 0
            t[5] = 0
            t[6] = 0
            t[7] = 0
            t[8] = 0
            t[9] = 0
            t[10] = 0
            t[11] = 0
            t[12] = 0
            t[13] = 0
            t[14] = 0
            t[15] = 0
            t[16] = 0
            t[17] = 0
            t[18] = 0
            t[19] = 0
            t[20] = 0
            t[21] = 0
          wt.append([i, organism_zip_file_name, no_of_proteins, t[0], t[1], t[2], t[3], t[4], t[5], t[6], t[7], t[8], t[9], t[10], t[11], t[12], t[13], t[14], t[15], t[16], t[17], t[18], t[19], t[20], t[21]])
    wf.save(xl_file_name_to_save)
    print(xl_file_name_to_save+" is saved succesfully!")

if __name__ == '__main__':
    main()

