# PACKAGES #

import os
import gzip
import subprocess
import pandas as pd
import zipfile

########### VARIABLES ##############

# Path where genomes have been downloaded
source_path = 'your_path'
# Path where the genomes and files will be stored after extraction
genomes_path = 'your_path'
# Path where retrozymes will be stored
retrozymes_path = 'your_path'
# Path where retrozymes sequences will be stored
retrozymes_seqs_path = 'your_path'

# Path to the ribozyme model 
ribozyme_model = 'your_path' 

# Max retrozyme length
rtzm_length = 2000

####################################


#### DECOMPRESSION AND CMSEARCH ####

# FUNCTION: Decompress ZIP files and handle errors
def decompress_zip(zipfile, destination_folder):
    try:
        # Verify if the ZIP file exists
        if not os.path.exists(zipfile):
            return None
        # Verify if the destination folder exists, if not, create it
        if not os.path.exists(destination_folder):
            os.makedirs(destination_folder)
        
        # Decompress the ZIP file
        with zipfile.ZipFile(zipfile, 'r') as zip_ref:
            zip_ref.extractall(destination_folder)
    except zipfile.BadZipFile:
        print("File " + zipfile + " is not a valid ZIP file.")
        if os.path.exists(destination_folder) and not os.listdir(destination_folder):
            os.rmdir(destination_folder)
        raise
    except Exception as e:
        print("Error while decompressing " + zipfile + ": " + str(e))
        if os.path.exists(destination_folder) and not os.listdir(destination_folder):
            os.rmdir(destination_folder)
        raise


# Retrieve assembly names and store them in a list
folder_list = os.listdir(source_path)
extension = '.zip'
folder_list_filtered = []
for element in folder_list:
    if element.endswith(extension):
        folder_list_filtered.append(element[:-element.index(extension)])  # Remove the .zip extension
genomes_list = folder_list_filtered.sort() # Sort them from lower to higher

# FUNCTION: Retrieve the path of the .fna file for a given genome
def fna_path(genome):
    fna_list = os.listdir(genomes_path + genome + '/ncbi_dataset/data/' + genome)
    return os.path.abspath(genomes_path + genome + '/ncbi_dataset/data/' + genome + '/' + fna_list[0])
    
# FUNCTION: Remove the .fna file and its directory to avoid storage issues   
def delete_fna(directory):
    for work_file in os.listdir(directory):
        path = os.path.join(directory, work_file)
        if os.path.isfile(path) or os.path.islink(path):
            os.unlink(path)
        elif os.path.isdir(path):
            os.rmdir(path)

# FUNCTION: Extract the description from the FASTA header           
def retrieve_description(fasta_path):
    with open(fasta_path, 'r') as f:
        for line in f:
            if line.startswith('>'):
                parts = line.strip().split(maxsplit=1) # Split in two parts: header and description
                if len(parts) > 1:  # Verify if there is a description
                    return parts[1]  # Return the description part
                else:
                    return None
    return None

# Create the directory for genomes and execute the CMSEARCH command for ribozyme search
for element in genomes_list:
    if os.path.exists(genomes_path + element):
        continue
    else:
        try:
            decompress_zip(source_path + element + '_dataset.zip', genomes_path + element)
        except Exception as e:
            print("Error while decompressing " + element + ": " + str(e))
            continue
    
        fa_description = retrieve_description(fna_path(element))
    
        # Execute the CMSEARCH command
        command = (
            'cmsearch '
            '-o ' + os.path.join(genomes_path, element, element + '_output.txt ') + # Output file
            '-A ' + os.path.join(genomes_path, element, element + '_alignments.txt ') + # Alignments file
            '--tblout ' + os.path.join(genomes_path, element, element + '_hits.txt ') + # Hits file
            ribozyme_model + ' ' +
            fna_path(element)
            )
        try:
            subprocess.run(command, shell=True, check=True)
        except subprocess.CalledProcessError:
            continue
    
        #### RETROZYMES ####
        
        # File processing: remove comments and create a DataFrame
        with open(os.path.join(genomes_path, element, element + '_hits.txt'), 'r') as work_file:
            lines = work_file.readlines()
            lines_without_comments = [line for line in lines if not line.strip().startswith('#')]
        with open(os.path.join(genomes_path, element, element + '_processed.txt'), 'w') as work_file:
            work_file.writelines(lines_without_comments)
        with open(os.path.join(genomes_path, element, element + '_processed.txt'), 'r') as work_file:
            data = [line.strip().split(None, 17) for line in work_file]
        hits_df = pd.DataFrame(data, columns=[
            'target_name', 'accession', 'query_name', 'accesion', 'mdl',
            'mdl_from', 'mdl_to', 'seq_from', 'seq_to', 'strand', 'trunc',
            'pass', 'gc', 'bias', 'score', 'E-value', 'inc', 'description'
        ])
        hits_df['description'] = hits_df['description'].str.replace(' ', '$', regex=False)
        hits_df['seq_from'] = pd.to_numeric(hits_df['seq_from'])
        hits_df['seq_to'] = pd.to_numeric(hits_df['seq_to'])
        hits_df = hits_df.sort_values(by=['target_name', 'seq_from'])
        hits_df.to_csv(os.path.join(genomes_path, element, element + '_processed.txt'), sep='\t', index=False)
    
        # Adjust seq_from and seq_to based on strand to take circRNA into account
        hits_df.loc[hits_df['strand'] == '+', 'seq_from'] += 5
        hits_df.loc[hits_df['strand'] == '-', 'seq_from'] -= 4
    
        # ---- Retrozymes (- strand) ----
            
            # Retrieve significant retrozymes in the negative strand
        hits_sign = hits_df[hits_df.inc == '!']
        hits_sign_str = hits_sign[hits_sign.strand == '-'].copy()
    
        # Sort by target_name and seq_from to ensure correct comparison
        hits_sign_str = hits_sign_str.sort_values(by=["target_name", "seq_from"]).reset_index(drop=True)
    
        # Initialize the 'distance' column with default values
        hits_sign_str['distance'] = 0
    
        # Calculate the distance for consecutive rows within the same target_name
        for i in range(1, len(hits_sign_str)):
            if hits_sign_str.loc[i, 'target_name'] == hits_sign_str.loc[i - 1, 'target_name']:
                hits_sign_str.at[i, 'distance'] = hits_sign_str.loc[i, 'seq_from'] - hits_sign_str.loc[i - 1, 'seq_from']
    
        # Filter rows where the distance is greater than 0 and less than or equal to the length threshold
        df_rows = []
        for i in range(len(hits_sign_str) - 1):
            if (0 < hits_sign_str.loc[i + 1, 'distance'] <= rtzm_length and hits_sign_str.loc[i, 'target_name'] == hits_sign_str.loc[i + 1, 'target_name']): 
                df_rows.append(hits_sign_str.iloc[i])
                df_rows.append(hits_sign_str.iloc[i + 1])
        # If the distance is greater than 0 and less than or equal to the length threshold, 
        # and both rows are included in the same contig, adds them to the list.
    
        # Create a DataFrame with the selected rows
        if df_rows:
            retrozymes = pd.DataFrame(df_rows).drop_duplicates().reset_index(drop=True)
            retrozymes.to_csv(retrozymes_path + element + '_retrozymes_neg.txt', sep=' ', index=False)
            
        # ---- Retrozymes (+ strand) ----
            
        # Retrieve significant retrozymes in the negative strand
        hits_sign = hits_df[hits_df.inc == '!']
        hits_sign_str = hits_sign[hits_sign.strand == '+'].copy()
    
        # Sort by target_name and seq_from to ensure correct comparison
        hits_sign_str = hits_sign_str.sort_values(by=["target_name", "seq_from"]).reset_index(drop=True)
    
        # Initialize the 'distance' column with default values
        hits_sign_str['distance'] = 0
    
        # Calculate the distance for consecutive rows within the same target_name
        for i in range(1, len(hits_sign_str)):
            if hits_sign_str.loc[i, 'target_name'] == hits_sign_str.loc[i - 1, 'target_name']:
                hits_sign_str.at[i, 'distance'] = hits_sign_str.loc[i, 'seq_from'] - hits_sign_str.loc[i - 1, 'seq_from']
    
        # Filter rows where the distance is greater than 0 and less than or equal to the length threshold
        df_rows = []
        for i in range(len(hits_sign_str) - 1):
            if (0 < hits_sign_str.loc[i + 1, 'distance'] <= rtzm_length and
                hits_sign_str.loc[i, 'target_name'] == hits_sign_str.loc[i + 1, 'target_name']):
                df_rows.append(hits_sign_str.iloc[i])
                df_rows.append(hits_sign_str.iloc[i + 1])
    
        # Create a DataFrame with the selected rows
        if df_rows:
            retrozymes = pd.DataFrame(df_rows).drop_duplicates().reset_index(drop=True)
            retrozymes.to_csv(retrozymes_path + element + '_retrozymes_pos.txt', sep=' ', index=False)
    
        # ---- Join files ----
        
        neg_file = os.path.join(retrozymes_path, element + '_retrozymes_neg.txt')
        pos_file = os.path.join(retrozymes_path, element + '_retrozymes_pos.txt')
        combined_file = os.path.join(retrozymes_path, element + '_retrozymes.txt')
    
        if os.path.exists(neg_file) and os.path.exists(pos_file):
            with open(combined_file, 'w') as cat:
                with open(pos_file, 'r') as pos:
                    cat.write(pos.read())
                with open(neg_file, 'r') as neg:
                    neg_lines = neg.readlines()[1:]
                    cat.writelines(neg_lines)
                os.remove(pos_file)
                os.remove(neg_file)
        elif os.path.exists(neg_file):
            os.rename(neg_file, combined_file)
        elif os.path.exists(pos_file):
            os.rename(pos_file, combined_file)
    
    #### BEDTOOLS ####
    
        if os.path.exists(retrozymes_path + element + '_retrozymes.txt'):
            with open(retrozymes_path + element + '_retrozymes.txt', 'r') as work_file:
                data = [line.strip().split(None, 18) for line in work_file]
        
            header = data[0]
            data = data[1:] 
            df = pd.DataFrame(data, columns=header)
                
            # Convert columns to numeric and clean up the description
            df['seq_from'] = pd.to_numeric(df['seq_from'], errors='coerce')
            df['seq_to'] = pd.to_numeric(df['seq_to'], errors='coerce')
            df['distance'] = pd.to_numeric(df['distance'], errors='coerce')
            df['description'] = df['description'].str.replace('$', ' ', regex=False)
                
            # Replace NaN values in 'distance' and distances greater than threshold with -1 
            df['distance'] = df['distance'].fillna(-1)
            df.loc[df['distance'] > 2000, 'distance'] = -1
                
            # Extract coordinates for bedtools input
            bedtools_input_list = []
            last_seq_from = None
            last_target_name = None
            last_strand = None
                
            for i in range(len(df) - 1):  # To avoid index out of range
                target_name = df.loc[i, 'target_name']
                seq_from = df.loc[i, 'seq_from']
                seq_from_next = df.loc[i + 1, 'seq_from']
                distance = df.loc[i, 'distance']
                distance_next = df.loc[i + 1, 'distance']
                strand = df.loc[i, 'strand']
                description = df.loc[i, 'description']        
                
                if distance_next <= 0: #
                    continue
                else: 
                    bedtools_input_list.append([target_name + ' ' + description, int(seq_from) - 1, int(seq_from_next) - 1, '.', '0', strand])

             # Save bedtools input to a DataFrame and write to a BED file
            bedtools_input = pd.DataFrame(bedtools_input_list, columns=['target_name', 'seq_from_1', 'seq_from_2', 'name', 'orf', 'strand'])
            bedtools_input.to_csv(retrozymes_seqs_path + element + '_input.bed', sep='\t', index=False, header=None)       
        
            if os.path.exists(retrozymes_seqs_path + element + '_input.bed'):
                          input_fa = fna_path(element)
                          input_bed = retrozymes_seqs_path + element + '_input.bed'
                          output_fa = retrozymes_seqs_path + element + '_output.fa'
                          command = ['bedtools', 'getfasta', '-fullHeader', '-s', '-fi', input_fa, '-bed', input_bed, '-fo', output_fa]
                          
                          command = (
                          'bedtools getfasta -fullHeader -s -fi ' +
                          input_fa +
                          ' -bed ' +
                          input_bed +
                          ' -fo ' +
                          output_fa
                          )
                          result = subprocess.run(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                          print("STDOUT:", result.stdout.decode())
                          print("STDERR:", result.stderr.decode())

                          
                          delete_fna(os.path.join(genomes_path, element, 'ncbi_dataset/data', element))
        else:
            delete_fna(os.path.join(genomes_path, element, 'ncbi_dataset/data', element))
            continue
            
    
            
    
