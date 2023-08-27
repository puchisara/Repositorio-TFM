#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import argparse
import pandas as pd
import logging
import os
import numpy as np

def parse_vcf():
    '''
    Download the VCF file containing the problematic sites of SARS-CoV-2.
    '''
    url = "https://raw.githubusercontent.com/W-L/ProblematicSites_SARS-CoV2/master/archived_vcf/problematic_sites_sarsCov2.2021-10-27-18%3A39.vcf"
    filename = "problematic_sites_sarsCov2.vcf"

    # Check if the file already exists
    if os.path.exists(filename):
        logging.info(f"File '{filename}' already exists. Skipping download.")
    else:
        # Download the file
        vcf = pd.read_csv(
            url,
            delim_whitespace=True,
            comment="#",
            names=("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO")
        )
        vcf.to_csv(filename, index=False)
        logging.info(f"File '{filename}' downloaded.")

    # Read the file
    vcf = pd.read_csv(filename)
    positions = tuple(vcf.loc[vcf.FILTER.isin(["mask"]), "POS"])
    return positions

def parse_file(file_path):
    '''
    Parse the TSV file containing the metadata of the patients.
    '''
    df = pd.read_csv(file_path, sep='\t') # Read the TSV file

    # Sort each group by 'Day' in descending order and then convert to dictionary
    data_dict = df.groupby('PATIENT_ID').apply(lambda g: g.sort_values('Day', ascending=True).drop('PATIENT_ID', axis=1).to_dict(orient='records')).to_dict()

    return data_dict

def parse_arguments():
    '''
    Parse the command-line arguments.
    '''
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('-i', '--input', required=True, type=str, help='Path to input TSV file.')
    return parser.parse_args()

def create_consolidated_data_frame(TABLE_DAYS):
    '''
    Create a Pandas DataFrame
    '''
    consolidated_data = pd.DataFrame(columns=['PATIENT', 'POS', 'REF', 'ALT'])
    for day in TABLE_DAYS:
        consolidated_data['ALTFREQ_'+str(day)] = None
        consolidated_data['COV_'+str(day)] = None
        consolidated_data['DEPTH_'+str(day)] = None
        consolidated_data['TOTAL_DP_'+str(day)] = None

    return consolidated_data

def process_single_file(file_path, label, day_single, positions):
    '''
    Process a single sample TSV file and return a Pandas DataFrame.
    '''
   df = pd.read_csv(file_path, sep='\t')

    depth_filter = df['TOTAL_DP'] >= 30
    freq_filter  = df['ALT_FREQ'] >= 0.05
    pass_filter = df['PASS'] == True
    num_reads = df['ALT_DP'] >= 2  
    num_reads_rev = df['ALT_RV'] >= 2
    indel_filter = ~df['ALT'].str.contains('[+-]', regex=True)
    num_reads_filter = indel_filter & (num_reads & num_reads_rev)

    df = df[freq_filter & depth_filter & pass_filter & (~indel_filter | num_reads_filter) & ~df['POS'].isin(positions)].copy()

    df = df[['POS', 'REF', 'ALT', 'ALT_FREQ', 'TOTAL_DP']].rename(columns={'ALT_FREQ': f'ALT_FREQ_{day_single}','TOTAL_DP': f'TOTAL_DP_{day_single}'})    df[f'COV_{day_single}'] = label
    
    return df

def merge_data_frames(consolidated_data, df):
    '''
    Merge the consolidated data frame with the data frame of a single sample.
    '''
    if consolidated_data.empty:
        consolidated_data = df
    else:
        consolidated_data = pd.merge(consolidated_data, df, on=['PATIENT', 'POS', 'REF', 'ALT'], how='outer')
    
    return consolidated_data

def fill_and_process_depth(consolidated_data, days, patient_cov):
    '''
    Fill the missing values in the consolidated data frame and process the depth files.
    '''
    consolidated_data = consolidated_data.fillna(0)
    for day_single, path_cov in zip(days, patient_cov):
        depth_df = pd.read_csv(path_cov, sep='\t', header=None, names=['REF', 'POS', 'DEPTH'])
        depth_df = depth_df[['POS', 'DEPTH']].rename(columns={'DEPTH': f'DEPTH_{day_single}'})
        consolidated_data = pd.merge(consolidated_data, depth_df, on=['POS'], how='inner')
        consolidated_data.loc[consolidated_data[f'DEPTH_{day_single}'] < 30, f'ALT_FREQ_{day_single}'] = np.nan

    return consolidated_data

def calculate_mutations_and_reversions(consolidated_data, days):
    '''
    Calculate the number of mutations and reversions between consecutive days.
    '''
    day_pairs = list(zip(days[:-1], days[1:]))
    for day1, day2 in day_pairs:
        consolidated_data[f'MUTATIONS_{day1}-{day2}'] = 0
        consolidated_data[f'REVERSIONS_{day1}-{day2}'] = 0

    for i, row in consolidated_data.iterrows():
        for day1, day2 in day_pairs:
            freq1 = row[f'ALT_FREQ_{day1}']
            freq2 = row[f'ALT_FREQ_{day2}']

            if freq1 < 0.05 and freq2 >= 0.05:
                consolidated_data.at[i, f'MUTATIONS_{day1}-{day2}'] += 1
            elif freq1 >= 0.05 and freq2 < 0.05:
                consolidated_data.at[i, f'REVERSIONS_{day1}-{day2}'] += 1

    return consolidated_data

def process_patient(patient_data:list, labels:list, positions:tuple, patient_id:str, days: list, TABLE_DAYS:list, patient_cov:list):
    '''
    Process the data of a single patient.
    '''
    consolidated_data = create_consolidated_data_frame(TABLE_DAYS)

    for file_path, label, day_single, path_cov in zip(patient_data, labels, days, patient_cov):
        df = process_single_file(file_path, label, day_single, positions)
        df['PATIENT'] = patient_id
        consolidated_data = merge_data_frames(consolidated_data, df)
        
    for label, day_single in zip(labels, days):
        consolidated_data[f'COV_{day_single}'] = label

    consolidated_data = fill_and_process_depth(consolidated_data, days, patient_cov)

    consolidated_data = calculate_mutations_and_reversions(consolidated_data, days)

    return consolidated_data

def main():
    #python3 get_snp.py -i f.tsv
    args = parse_arguments()
    positions = parse_vcf()
    data = parse_file(args.input)
    TABLE_DAYS = pd.read_csv(args.input, sep='\t')['Day'].unique().tolist()

    all_data = []  # list to store all dataframes

    for patient_id, records in data.items():
        patient_data = [os.path.join(record["PATH_TSV"], record["COV_ID"] + ".tsv") for record in records] # list of paths to all tsv files
        patient_cov = [os.path.join(record["PATH_DEP"], record["COV_ID"] + ".depth") for record in records] # list of paths to all depth files
        labels = [record["COV_ID"] for record in records] # list of labels
        day = [record["Day"] for record in records] # list of days
        consolidated_data = process_patient(patient_data, labels, positions, patient_id, day, TABLE_DAYS, patient_cov) # process the data of a single patient
        all_data.append(consolidated_data)  # append the dataframe to the list

    # concatenate all dataframes
    all_data_df = pd.concat(all_data, ignore_index=True)

    # save the merged dataframe to a tsv file
    all_data_df.to_csv(f'{args.input}all_consolidated_data.tsv', sep='\t', index=False)

if __name__ == '__main__':
    main()
