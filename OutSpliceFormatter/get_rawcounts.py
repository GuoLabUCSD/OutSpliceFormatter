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

#Function to extract the number of unique reads from every Log.final.out file in Tumor and Normal STAR Directories

def combine_UMR(star_results_tumors, star_results_normals):
    combined_df = pd.DataFrame(columns = ['Descriptions'])
    
    for file in os.listdir(star_results_tumors):
        star_results_df = pd.DataFrame()
        if file.endswith("Log.final.out"):
            path = os.path.join(star_results_tumors, file)
            results = pd.read_csv(path, sep = '\t')
            star_results_df = star_results_df.append(results)
            star_results_df.columns = ['Descriptions', '{}'.format(file)]
            star_results_df = star_results_df.filter(items = [6], axis=0)
            filename = file.split('Log')[0]
            #filename = filename.ljust(10, 'X')
            star_results_df.rename(columns = {'{}'.format(file):'{}'.format(filename)}, inplace = True)
            combined_df = combined_df.merge(star_results_df, on = 'Descriptions', how = 'outer')
            
    for file in os.listdir(star_results_normals):
        star_results_df = pd.DataFrame()
        if file.endswith("Log.final.out"):
            path = os.path.join(star_results_normals, file)
            results = pd.read_csv(path, sep = '\t')
            star_results_df = star_results_df.append(results)
            star_results_df.columns = ['Descriptions', '{}'.format(file)]
            star_results_df = star_results_df.filter(items = [6], axis=0)
            filename = file.split('Log')[0]
            #filename = filename.ljust(10, 'X')
            star_results_df.rename(columns = {'{}'.format(file):'{}'.format(filename)}, inplace = True)
            combined_df = combined_df.merge(star_results_df, on = 'Descriptions', how = 'outer')
    
    combined_df = combined_df.drop(['Descriptions'], axis = 1)
    combined_df = combined_df.reindex(sorted(combined_df.columns), axis=1)
    combined_df.insert(0, 'Description', 'Number of Uniquely Mapped Reads')
    return combined_df

A = combine_UMR(args.star_results_tumors, args.star_results_normals)
A.to_csv('{}/{}_rawcounts.txt'.format(args.output_directory, args.output_prefix), sep = '\t', index = False)

