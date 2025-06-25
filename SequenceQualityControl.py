# PACKAGES

import re
import os

#### Path definitions ####

QC_rtzms = "your_path" # Path to the directory containing the retrozyme files
QC_characters = "your_path" #
QC_Ns = "your_path"         # Paths where the processed 
QC_length = "your_path"     # files will be saved
QC_LTRs = "your_path"       #

#### Replaces invalid characters in FASTA entries for tree building ####

def clean_text(input_file, output_file):
    with open(input_file, "r") as infile, open(output_file, "w") as outfile:
        for line in infile:
            line = re.sub(r"\(\+\)", " + ", line)  # Replaces '(+)' with ' + '
            line = re.sub(r"\(-\)", " - ", line)  # Replaces '(-)' with ' - '
            line = line.replace(":", " ")  # Replaces ':' with ' '
            line = line.replace(",", " ")  # Replaces ',' with ' '
            outfile.write(line)

file_list = os.listdir(QC_rtzms)

for element in file_list:
    input_fasta = QC_rtzms + element
    output_fasta = QC_characters + element
    clean_text(input_fasta, output_fasta)



#### Delete entries with unknown nucleotides in their sequences (Ns) ####

def filter_fasta_Ns(input_fasta, output_fasta):
    with open(input_fasta, "r") as infile, open(output_fasta, "w") as outfile:
        keep_entry = False
        header = ""
        sequence = ""
        
        for line in infile:
            if line.startswith(">"):
                if keep_entry and sequence:
                    outfile.write(header + sequence + "\n")
                header = line
                sequence = ""
                keep_entry = True
            else:
                sequence += line.strip()
                if "N" in line:
                    keep_entry = False
        
        if keep_entry and sequence:
            outfile.write(header + sequence + "\n")

file_list = os.listdir(QC_characters)

for element in file_list:
    input_fasta = QC_characters + element
    output_fasta = QC_Ns + element
    filter_fasta_Ns(input_fasta, output_fasta)



#### Delete retrozymes smaller or larger than the given values (inactive or polymers) ####

filter_smaller_than = 500
filter_bigger_than = 1150


def filter_fasta_polymers(input_file, output_file, min_length=filter_smaller_than, max_length=filter_bigger_than):
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        fasta_entries = infile.read().strip().split('>')
        
        for entry in fasta_entries:
            if not entry:
                continue
            
            lines = entry.split('\n')
            header = lines[0]
            sequence = ''.join(lines[1:])
            
            match = re.search(r"(\d+)-(\d+)", header)
            if match:
                start = int(match.group(1))
                end = int(match.group(2))
                length = end - start
                
                if length > min_length and length < max_length:
                    outfile.write('>' + header + '\n' + sequence + '\n')


file_list = os.listdir(QC_Ns)

for element in file_list:
    input_fasta = QC_Ns + element
    output_fasta = QC_length + element
    filter_fasta_polymers(input_fasta, output_fasta)
    
    
    
#### Extracts the sequence that includes the 2 LTRs of the retrozymes (5'-260 + 260-3') ####

segment_length = 260  # Length of the segments to extract from each end of the sequence

def extract_fasta_ltrs(input_file, output_file, segment_length):
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        fasta_entries = infile.read().strip().split('>')
        
        for entry in fasta_entries:
            if not entry:
                continue
            
            lines = entry.split('\n')
            header = lines[0]
            sequence = ''.join(lines[1:])
            
            if len(sequence) >= segment_length * 2:
                extracted_sequence = sequence[:segment_length] + sequence[-segment_length:]
            else:
                extracted_sequence = sequence  # If the sequence is shorter than 520, keep it as is
            
            outfile.write('>' + header + '\n' + extracted_sequence + '\n')

file_list = os.listdir(QC_length)

for element in file_list:
    input_fasta = QC_length + element
    output_fasta = QC_LTRs + element
    extract_fasta_ltrs(input_fasta, output_fasta)    