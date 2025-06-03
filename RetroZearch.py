# -*- coding: utf-8 -*-

import os
import gzip
import subprocess
import pandas as pd
import zipfile

###########VARIABLES##############

genomes_path = os.getcwd() + '/3er_analisis/1.genomes/'
rivero_path = '/media/disk/home/rivero/genomes/spermatophyta_NCBI_cntg_scff/chrom_fulll_Magnoliophyta/'

longitud_rtzm = 2000

##################################


#### DESCOMPRESION Y CMSEARCH ####

## Descomprimir archivos
def descomprimir_zip(archivo_zip, carpeta_destino):
    try:
        # Verifica si el archivo ZIP existe
        if not os.path.exists(archivo_zip):
            return None
        # Verifica si la carpeta de destino existe, si no, la crea
        if not os.path.exists(carpeta_destino):
            os.makedirs(carpeta_destino)
        
        # Descomprimir el archivo ZIP
        with zipfile.ZipFile(archivo_zip, 'r') as zip_ref:
            zip_ref.extractall(carpeta_destino)
    except zipfile.BadZipFile:
        print("El archivo " + archivo_zip + " no es un archivo ZIP valido.")
        if os.path.exists(carpeta_destino) and not os.listdir(carpeta_destino):
            os.rmdir(carpeta_destino)
        raise
    except Exception as e:
        print("Error al descomprimir " + archivo_zip + ": " + str(e))
        if os.path.exists(carpeta_destino) and not os.listdir(carpeta_destino):
            os.rmdir(carpeta_destino)
        raise


# Lista de genomas
root_path = genomes_path  # obtiene la ruta del script

folder_list = os.listdir(rivero_path)
extension = '.zip'
folder_list_filtered = []

for element in folder_list:
    if element.endswith(extension):
        folder_list_filtered.append(element[:-12])
folder_list_filtered.sort()

folder_list_2_analysis = os.listdir(os.getcwd() + '/preparacion/1.genomes2.0/')

genomes_list = [x for x in folder_list_filtered if x not in folder_list_2_analysis] # Quita los duplicados del 2o analisis


# Obtener ruta del .fna
def ruta_fna(genoma):
    fna_list = os.listdir(root_path + genoma + '/ncbi_dataset/data/' + genoma)
    return os.path.abspath(root_path + genoma + '/ncbi_dataset/data/' + genoma + '/' + fna_list[0])
    
# Eliminar el fasta    
def eliminar_fna(directorio):
    for archivo in os.listdir(directorio):
        ruta = os.path.join(directorio, archivo)
        if os.path.isfile(ruta) or os.path.islink(ruta):  # Archivos o enlaces simbolicos
            os.unlink(ruta)  # Elimina el archivo
        elif os.path.isdir(ruta):  # Subdirectorios
            os.rmdir(ruta)  # Elimina el subdirectorio vacio

# Extraer la descripción del fasta (especie)            
def extraer_descripcion(fasta_path):
    with open(fasta_path, 'r') as f:
        for line in f:
            if line.startswith('>'):  # Línea de encabezado
                partes = line.strip().split(maxsplit=1)  # Dividir solo en dos partes
                if len(partes) > 1:  # Verificar si hay descripción
                    return partes[1]  # Devolver la descripción completa
                else:
                    return None  # Si no hay descripción
    return None  # Si el archivo está vacío o no tiene encabezados

# Creacion de carpeta y ejecucion del comando
for element in genomes_list:
    
    if os.path.exists(root_path + element):
        continue
    
    else:
        print(rivero_path + element + '.zip')
        try:
            descomprimir_zip(rivero_path + element + '_dataset.zip', root_path + element)
        except Exception as e:
            print("Error al descomprimir " + element + ": " + str(e))
            continue
    
        fa_description = extraer_descripcion(ruta_fna(element))
    
        # Ejecutar el comando
        command = (
            'cmsearch '
            '-o ' + os.path.join(root_path, element, element + '_output.txt ') +
            '-A ' + os.path.join(root_path, element, element + '_alignments.txt ') +
            '--tblout ' + os.path.join(root_path, element, element + '_hits.txt ') +
            'hhr_270322.cm ' + 
            ruta_fna(element)
            )
        try:
        # Ejecutar el comando
            subprocess.run(command, shell=True, check=True)
        except subprocess.CalledProcessError:
            continue
    
        #### RETROZIMAS ####
        
        with open(os.path.join(root_path, element, element + '_hits.txt'), 'r') as archivo:
            lines = archivo.readlines()
            lines_without_comments = [line for line in lines if not line.strip().startswith('#')]
        
        with open(os.path.join(root_path, element, element + '_processed.txt'), 'w') as archivo:
            archivo.writelines(lines_without_comments)
        
        with open(os.path.join(root_path, element, element + '_processed.txt'), 'r') as archivo:
            data = [line.strip().split(None, 17) for line in archivo]
        
        hits_df = pd.DataFrame(data, columns=[
            'target_name', 'accession', 'query_name', 'accesion', 'mdl',
            'mdl_from', 'mdl_to', 'seq_from', 'seq_to', 'strand', 'trunc',
            'pass', 'gc', 'bias', 'score', 'E-value', 'inc', 'description'
        ])
        
        hits_df['description'] = hits_df['description'].str.replace(' ', '$', regex=False)
        hits_df['seq_from'] = pd.to_numeric(hits_df['seq_from'])
        hits_df['seq_to'] = pd.to_numeric(hits_df['seq_to'])
        hits_df = hits_df.sort_values(by=['target_name', 'seq_from'])
        hits_df.to_csv(os.path.join(root_path, element, element + '_processed.txt'), sep='\t', index=False)
    
        # Desplazar 4pb el comienzo de la retrozima
        hits_df.loc[hits_df['strand'] == '+', 'seq_from'] += 5
        hits_df.loc[hits_df['strand'] == '-', 'seq_from'] -= 4
    
        # ---- Retrozimas (-) ----
            
            # Filtrar las retrozimas significativas con la cadena negativa
        hits_sign = hits_df[hits_df.inc == '!']
        hits_sign_str = hits_sign[hits_sign.strand == '-'].copy()
    
        # Ordenar por target_name y seq_from para asegurar la correcta comparación
        hits_sign_str = hits_sign_str.sort_values(by=["target_name", "seq_from"]).reset_index(drop=True)
    
        # Inicializar la columna 'distancia' con valores predeterminados
        hits_sign_str['distancia'] = 0
    
        # Calcular la distancia para filas consecutivas dentro del mismo target_name
        for i in range(1, len(hits_sign_str)):
            if hits_sign_str.loc[i, 'target_name'] == hits_sign_str.loc[i - 1, 'target_name']:
                hits_sign_str.at[i, 'distancia'] = hits_sign_str.loc[i, 'seq_from'] - hits_sign_str.loc[i - 1, 'seq_from']
    
        # Filtrar filas donde la distancia sea mayor que 0 y menor o igual a 2000
        filas = []
        for i in range(len(hits_sign_str) - 1):
            if (0 < hits_sign_str.loc[i + 1, 'distancia'] <= longitud_rtzm and hits_sign_str.loc[i, 'target_name'] == hits_sign_str.loc[i + 1, 'target_name']): # Si la distancia entre dos ribozimas esta entre 0 y 2000, y forman parte del mismo contig, las añade al nuevo archivo
                filas.append(hits_sign_str.iloc[i])
                filas.append(hits_sign_str.iloc[i + 1])
    
        # Crear un DataFrame con las filas seleccionadas
        print(os.path.join(os.getcwd(), '/3er_analisis/2.retrozymes/', element + '_retrozymes_neg.txt'))
        if filas:
            retrozymes = pd.DataFrame(filas).drop_duplicates().reset_index(drop=True)
            retrozymes.to_csv(os.path.join(os.getcwd(), '3er_analisis/2.retrozymes', element + '_retrozymes_neg.txt'), sep=' ', index=False)
            
        # ---- Retrozimas (+) ----
            
        # Filtrar las retrozimas con la cadena negativa
        hits_sign = hits_df[hits_df.inc == '!']
        hits_sign_str = hits_sign[hits_sign.strand == '+'].copy()
    
        # Ordenar por target_name y seq_from para asegurar la correcta comparación
        hits_sign_str = hits_sign_str.sort_values(by=["target_name", "seq_from"]).reset_index(drop=True)
    
        # Inicializar la columna 'distancia' con valores predeterminados
        hits_sign_str['distancia'] = 0
    
        # Calcular la distancia para filas consecutivas dentro del mismo target_name
        for i in range(1, len(hits_sign_str)):
            if hits_sign_str.loc[i, 'target_name'] == hits_sign_str.loc[i - 1, 'target_name']:
                hits_sign_str.at[i, 'distancia'] = hits_sign_str.loc[i, 'seq_from'] - hits_sign_str.loc[i - 1, 'seq_from']
    
        # Filtrar filas donde la distancia sea mayor que 0 y menor o igual a 2000
        filas = []
        for i in range(len(hits_sign_str) - 1):
            if (0 < hits_sign_str.loc[i + 1, 'distancia'] <= longitud_rtzm and
                hits_sign_str.loc[i, 'target_name'] == hits_sign_str.loc[i + 1, 'target_name']):
                filas.append(hits_sign_str.iloc[i])
                filas.append(hits_sign_str.iloc[i + 1])
    
        # Crear un DataFrame con las filas seleccionadas
        if filas:
            retrozymes = pd.DataFrame(filas).drop_duplicates().reset_index(drop=True)
            retrozymes.to_csv(os.path.join(os.getcwd(), '3er_analisis/2.retrozymes', element + '_retrozymes_pos.txt'), sep=' ', index=False)
    
        # ---- Concatenar archivos ----
        rz_path = os.getcwd() + '/3er_analisis/2.retrozymes/'
        
        neg_file = os.path.join(rz_path, element + '_retrozymes_neg.txt')
        pos_file = os.path.join(rz_path, element + '_retrozymes_pos.txt')
        combined_file = os.path.join(rz_path, element + '_retrozymes.txt')
    
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
    
    #### SUBSEQ ####
    
        if os.path.exists(os.getcwd() + '/3er_analisis/2.retrozymes/' + element + '_retrozymes.txt'):
            with open(os.getcwd() + '/3er_analisis/2.retrozymes/' + element + '_retrozymes.txt', 'r') as archivo:
                data = [line.strip().split(None, 18) for line in archivo]
        
            encabezado = data[0]
            data = data[1:]
            df = pd.DataFrame(data, columns=encabezado)
                
            # Convertir columnas numéricas
            df['seq_from'] = pd.to_numeric(df['seq_from'], errors='coerce')
            df['seq_to'] = pd.to_numeric(df['seq_to'], errors='coerce')
            df['distancia'] = pd.to_numeric(df['distancia'], errors='coerce')
            df['description'] = df['description'].str.replace('$', ' ', regex=False)
                
            # Reemplazar valores anómalos o NaN
            df['distancia'] = df['distancia'].fillna(-1)
            df.loc[df['distancia'] > 2000, 'distancia'] = -1
                
            # Extraer las secuencias para subseq
            subseq_input_list = []
            last_seq_from = None
            last_target_name = None
            last_strand = None
                
            for i in range(len(df) - 1):  # Esto evita el indice fuera de rango
                target_name = df.loc[i, 'target_name']
                seq_from = df.loc[i, 'seq_from']
                seq_from_next = df.loc[i + 1, 'seq_from']
                distancia = df.loc[i, 'distancia']
                distancia_next = df.loc[i + 1, 'distancia']
                strand = df.loc[i, 'strand']
                description = df.loc[i, 'description']        
                
                if distancia_next <= 0:
                    continue
                else:
                    subseq_input_list.append([target_name + ' ' + description, int(seq_from) - 1, int(seq_from_next) - 1, '.', '0', strand])

              
            print(subseq_input_list)      
             # Guardar el archivo _input.txt
            subseq_input = pd.DataFrame(subseq_input_list, columns=['target_name', 'seq_from_1', 'seq_from_2', 'name', 'orf', 'strand'])
            subseq_input.to_csv(os.getcwd() + '/3er_analisis/3.1.retrozymes_seqs/' + element + '_input.bed', sep='\t', index=False, header=None)       
        
            if os.path.exists(os.getcwd() + '/3er_analisis/3.1.retrozymes_seqs/' + element + '_input.bed'):
                          print('existe directorio')
                          input_fa = ruta_fna(element)
                          input_bed = os.getcwd() + '/3er_analisis/3.1.retrozymes_seqs/' + element + '_input.bed'
                          output_fa = os.getcwd() + '/3er_analisis/3.1.retrozymes_seqs/' + element + '_output.fa'
                          print(output_fa)
                          command = ['bedtools', 'getfasta', '-fullHeader', '-s', '-fi', input_fa, '-bed', input_bed, '-fo', output_fa]
                          print(command)
                          
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

                          
                          eliminar_fna(os.path.join(os.getcwd(), '3er_analisis/1.genomes', element, 'ncbi_dataset/data', element))
        else:
            eliminar_fna(os.path.join(os.getcwd(), '3er_analisis/1.genomes', element, 'ncbi_dataset/data', element))
            continue
            
    
            
    