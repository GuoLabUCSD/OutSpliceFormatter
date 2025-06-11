import argparse
import pandas as pd
import os
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

parser = argparse.ArgumentParser()
args_path = parser.add_argument_group("Inputs:")

args_path.add_argument('-c', '--rsem_results_tumors', type=str, required=True, help='Path to Tumor Directory with genes.results files from RSEM')
args_path.add_argument('-o', '--output_directory', type=str, required=True, help='Location to store output files')
args_path.add_argument('-p', '--output_prefix', type=str, required=True, help='Prefix to name output files')

args = parser.parse_args()

#Function to get the expected counts from RSEM that will undergo normalization

def combine_rsem_rawcounts(rsem_results_tumors):
    combined_df = pd.DataFrame(columns = ['gene_id'])
    
    for file in os.listdir(rsem_results_tumors):
        rsem_results_df = pd.DataFrame()
        if file.endswith("genes.results"):
            path = os.path.join(rsem_results_tumors, file)
            results = pd.read_csv(path, sep = '\t')
            rsem_results_df = rsem_results_df.append(results)
            rsem_results_df = rsem_results_df.drop(['transcript_id(s)', 'length', 'effective_length', 'TPM', 'FPKM'], axis = 1)
            filename = file.split('.genes.')[0]
            #filename = filename.ljust(10, 'X')
            rsem_results_df.rename(columns = {'expected_count':'{}'.format(filename)}, inplace = True)
            combined_df = combined_df.merge(rsem_results_df, on = 'gene_id', how = 'outer')
    
    combined_df = combined_df.reindex(sorted(combined_df.columns), axis=1)
    raw_id_col = combined_df['gene_id']
    combined_df.drop(['gene_id'], axis=1, inplace = True)
    combined_df.insert(0, 'gene_id', raw_id_col)
    return combined_df


X = combine_rsem_rawcounts(args.rsem_results_tumors)
X.to_csv('{}/{}_expected_counts_TMP.txt'.format(args.output_directory, args.output_prefix), sep = '\t', index = False)

