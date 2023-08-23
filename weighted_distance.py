#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import pandas as pd
from Bio import SeqIO
import os
import logging
from typing import List
import numpy as np
import multiprocessing
from joblib import Parallel, delayed, effective_n_jobs
import argparse
from itertools import combinations


def read_and_filter_tsv(file, positions):
    df = pd.read_table(file)
    # Define filters
    depth_filter = df['TOTAL_DP'] >= 30  # Minimum read depth
    freq_filter = df['ALT_FREQ'] >= 0.05  # Minimum allele frequency
    pass_filter = df['PASS'] == True  # 'PASS' flag must be true
    num_reads = df['ALT_DP'] >= 2  # Minimum number of alternative reads
    num_reads_rev = df['ALT_RV'] >= 2  # Minimum number of reverse reads

    # Filter for rows without insertion/deletion
    indel_filter = ~df['ALT'].str.contains('[+-]', regex=True)

    # Define a composite filter
    num_reads_filter = indel_filter & (num_reads & num_reads_rev)

    # Apply filters to the data and append to the 'filt' list
    df = df[freq_filter & depth_filter & pass_filter & (~indel_filter | num_reads_filter) & ~df['POS'].isin(positions)]
    df['REGION'] = os.path.splitext(os.path.basename(file))[0]  # Add 'REGION' column
    return df

# Add reference to the analysis
def tsv_from_seq(tsv_reference, reference, reference_name):

    mask = parse_vcf()
    seq = tsv_reference
    reference = reference

    pos = []
    ALT = []

    for i in range(0, len(seq)):

        if i not in mask:

            pos.append(i + 1)
            ALT.append(reference[i])

    df = pd.DataFrame({"POS": pos, "ALT": ALT})
    df["ALT_FREQ"] = 1
    df["REGION"] = reference_name

    return df

def read_and_concatenate_tsvs(file_list, tsv_reference, reference, reference_name):
    df_list = []
    positions = parse_vcf()  # Get the mask positions

    for file in file_list:
        df = read_and_filter_tsv(file, positions)
        df_list.append(df)
    df_list.append(tsv_from_seq(tsv_reference, reference, reference_name))
    concatenated_df = pd.concat(df_list, ignore_index=True)
    return concatenated_df

def create_freq_dict(filepaths: List[str]) -> dict:
    alt_values = set()  # Utilizamos un conjunto para evitar duplicados

    # Itera a través de cada archivo y agrega los valores únicos de 'ALT' al conjunto
    for filepath in filepaths:
        df = pd.read_table(filepath)
        alt_values.update(df['ALT'].unique())

    # Crea el diccionario freq a partir del conjunto de valores de 'ALT', con un valor inicial de 0 para todas las claves
    freq = {alt_value: 0 for alt_value in alt_values}
    
    # Agrega la clave 'A' con valor cero si no está presente en el conjunto de valores únicos de 'ALT'
    if 'A' not in freq:
        freq['A'] = 0

    return freq

def parse_vcf():
    url = "https://raw.githubusercontent.com/W-L/ProblematicSites_SARS-CoV2/master/archived_vcf/problematic_sites_sarsCov2.2021-10-27-18%3A39.vcf"
    filename = "problematic_sites_sarsCov2.vcf"

    # Check if the file already exists
    if os.path.exists(filename):
        pass
    else:
        # Download the file
        vcf = pd.read_csv(
            url,
            delim_whitespace=True,
            comment="#",
            names=("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO")
        )
        vcf.to_csv(filename, index=False)

    # Read the file
    vcf = pd.read_csv(filename)
    positions = tuple(vcf.loc[vcf.FILTER.isin(["mask"]), "POS"])
    return positions

def get_pos_tup(df, pos, reference, freq):
    freq = freq.copy()  # Hace una copia de freq para no modificar el diccionario original
    alt_keys = sorted(list(freq.keys()))  # Crea una lista ordenada de las llaves

    # Si la posicion esta en el dataframe, obtiene las frecuencias alelicas
    if pos + 1 in df["POS"].values:
        df_ = df[df["POS"] == pos + 1]
        for base in alt_keys:
            if base in df_["ALT"].values:
                freq[base] = float(df_["ALT_FREQ"][df_["ALT"] == base].iloc[0])

    # Calcula la frecuencia de la base en la secuencia de referencia
    ref = 1 - sum(freq.values())
    freq[reference[pos]] += ref
    # Retorna una tupla de las frecuencias alelicas en el orden de alt_keys
    return tuple(freq[base] for base in alt_keys)

def heterozygosity(freqs):
    return 1 - sum([f ** 2 for f in freqs])

def calc_heterozygosities(df1, df2, pos, reference, freq):
    freqs1 = get_pos_tup(df1, pos, reference, freq)
    freqs2 = get_pos_tup(df2, pos, reference, freq)

    hs1 = heterozygosity(freqs1)
    hs2 = heterozygosity(freqs2)
    hs = (hs1 + hs2) / 2

    total_freqs = [(f1 + f2) / 2 for f1, f2 in zip(freqs1, freqs2)]
    ht = heterozygosity(total_freqs)
    return hs, ht

def calc_fst_weir_cockerham(hs, ht):
    return (ht - hs) / ht if ht != 0 else 0

def get_dif_n(df, COV1, COV2, mask, reference, freq, i):
    df1 = df[df["REGION"] == COV1]
    df2 = df[df["REGION"] == COV2]
    result = calc_fst_weir_cockerham(*calc_heterozygosities(df1, df2, i, reference, freq))
    return result

def get_matrix(df, cov_list, reference, freq, directory=None):
    mask_positions = parse_vcf()

    def calculate_distance(sample1, sample2, directory):
        print(f"Calculating distances for samples: {sample1} and {sample2}")
        distances = []
        for i in range(len(reference)):
            if i + 1 not in mask_positions:
                distance1 = get_dif_n(df, sample1, sample2, mask_positions, reference, freq, i)
                distances.append((i+1, sample1, distance1))  # Agregar la posición, el nombre de la muestra y la distancia
        with open(os.path.join(directory, f"{sample1}_{sample2}_distances.txt"), "w") as file:
            for distance in distances:
                file.write(f"Position: {distance[0]}, Sample1: {distance[1]}, Distance: {distance[2]}\n")
            print(f"Distances for samples {sample1} and {sample2} saved in {sample1}_{sample2}_distances.txt")

    num_samples = len(cov_list)
    directory = directory or "/home/puchi/pos_ponde/"

    print(f"Generating distance files for all sample pairs...")
    sample_pairs = combinations(cov_list, 2)  # Obtener todas las combinaciones únicas de pares de muestras
    for sample1, sample2 in sample_pairs:
        calculate_distance(sample1, sample2, directory)
    

    for i, sample in enumerate(cov_list):
        distance_matrix[sample] = [result[i] for result in results]

    return pd.DataFrame(distance_matrix)

def read_input_file(input_file):
    # Lee el archivo .tsv en un DataFrame de pandas
    input_df = pd.read_table(input_file, sep='\t')
    # Combina las columnas 'path' e 'id' para crear la lista de archivos .tsv
    file_list = (input_df['path'] + input_df['id'] + ".tsv").tolist()
    return file_list

def parse_arguments():
    parser = argparse.ArgumentParser(description='Process TSV files and generate a distance matrix.')
    parser.add_argument('-i', '--input', required=True, type=str, help='Path to input TSV file.')
    parser.add_argument('-rt', '--reference_tsv', required=True, type=str, help='Path to reference FASTA file for tsv.')
    parser.add_argument('-r', '--reference', required=True, type=str, help='Path to reference FASTA file(Could be the same as --reference_tsv).')
    parser.add_argument('-o', '--output', required=True, type=str, help='Path to output TXT file.')
    return parser.parse_args()

def main():
    args = parse_arguments()
    input_file = args.input
    reference = str(next(SeqIO.parse(args.reference_tsv, "fasta")).seq)
    outgroup = str(next(SeqIO.parse(args.reference, "fasta")).seq)
    outgroup_name = str(next(SeqIO.parse(args.reference, "fasta")).id)
    file_list = read_input_file(input_file)
    df = read_and_concatenate_tsvs(file_list, reference, outgroup, outgroup_name)
    cov_list = df["REGION"].unique().tolist()
    freq = create_freq_dict(file_list)
    df = get_matrix(df, cov_list, reference, freq)
    np.savetxt(args.output, df.values, fmt='%.6f', delimiter='\t', header='\t'.join(df.columns))
    df = df.astype(float)
    df.columns = df.columns.astype(int)
    df.to_csv(args.output + ".csv")

if __name__ == "__main__":
    main()
