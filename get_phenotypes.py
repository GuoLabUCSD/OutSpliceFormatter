import argparse
import pandas as pd
import os
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

parser = argparse.ArgumentParser()
args_path = parser.add_argument_group("Inputs:")

args_path.add_argument('-a', '--star_results_tumors', type=str, required=True, help='Path to Tumor Directory with Log.final.out files from STAR')
args_path.add_argument('-b', '--star_results_normals', type=str, required=True, help='Path to Normal Directory with Log.final.out files from STAR')
args_path.add_argument('-o', '--output_directory', type=str, required=True, help='Location to store output files')
args_path.add_argument('-p', '--output_prefix', type=str, required=True, help='Prefix to name output files')

args = parser.parse_args()

#Function to get Phenotype Matrix based on given directories

def get_phenotypes(tumor_dir, normal_dir):
    pheno_matrix = pd.DataFrame(columns = ['Sample', 'Pheno'])
    
    index = 0
    for file in os.listdir(tumor_dir):
        if file.endswith("SJ.out.tab"):
            filename = file.split('SJ.out.')[0]
            pheno_matrix.loc[index] = [filename, 'T']
            index = index + 1
            
    for file in os.listdir(normal_dir):
        if file.endswith("SJ.out.tab"):
            filename = file.split('SJ.out.')[0]
            pheno_matrix.loc[index] = [filename, 'F']
            index = index + 1
            
    pheno_matrix = pheno_matrix.sort_values('Sample')
    return pheno_matrix

P = get_phenotypes(args.star_results_tumors, args.star_results_normals)
P.to_csv('{}/{}_phenos.txt'.format(args.output_directory, args.output_prefix), sep = '\t', index = False)
P

