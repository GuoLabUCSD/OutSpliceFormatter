import argparse
import pandas as pd
import os
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

parser = argparse.ArgumentParser()
args_path = parser.add_argument_group("Inputs:")

args_path.add_argument('-o', '--output_directory', type=str, required=True, help='Location to store output files')
args_path.add_argument('-p', '--output_prefix', type=str, required=True, help='Prefix to name output files')

args = parser.parse_args()

#Function to Re-Format the normalized RSEM Data from TCGA V2 Pipeline

def format_normalized_rsem(path_to_normalized_rsem_table):
    normalized_results_df = pd.read_csv(path_to_normalized_rsem_table, sep = '\t')
    normalized_results_df.columns = normalized_results_df.columns.str.replace("[.]", "-")
    #normalized_results_df.loc[-1] = 'normalized_count'
    #normalized_results_df.index = normalized_results_df.index + 1
    #normalized_results_df = normalized_results_df.sort_index()

    return normalized_results_df

y = format_normalized_rsem('{}/{}_genes_normalized_TMP.txt'.format(args.output_directory, args.output_prefix))
y.to_csv('{}/{}_genes_normalized.txt'.format(args.output_directory, args.output_prefix), sep = '\t', index = False)

