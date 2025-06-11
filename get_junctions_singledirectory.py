import argparse
import pandas as pd
import os
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

parser = argparse.ArgumentParser()
args_path = parser.add_argument_group("Inputs:")

args_path.add_argument('-a', '--star_results_tumors', type=str, required=True, help='Path to Tumor Directory with Log.final.out files from STAR')
args_path.add_argument('-o', '--output_directory', type=str, required=True, help='Location to store output files')
args_path.add_argument('-p', '--output_prefix', type=str, required=True, help='Prefix to name output files')

args = parser.parse_args()

#Function to get Junction Counts from STAR SJ.out.tab

acceptable_values = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 
                     'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 
                     'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY']

def combine_junction_counts(star_results_tumors):
    combined_junction = pd.DataFrame(columns = ['junction'])

    for file in os.listdir(star_results_tumors):
        junction_counts_df = pd.DataFrame()
        formatted = False
        if file.endswith("SJ.out.tab"):
            path = os.path.join(star_results_tumors, file)
            junction_counts = pd.read_csv(path, sep = '\t', header = None, names = ['chromosome', 'intron_start', 'intron_end', 'strand', 'motif', 'is_annotated', 'unique_reads', 'multi_mapped_reads', 'max_overhang'])
            check_firstcol = junction_counts['chromosome'].to_list()
            for i in check_firstcol:
                i = str(i)
                if 'chr' in i:
                    formatted = True
                    break
            if formatted == False:
                junction_counts['chromosome'] = 'chr' + junction_counts['chromosome'].astype(str)
                check_firstcol = junction_counts['chromosome'].to_list()
            ind2lose = []
            for index, string in enumerate(check_firstcol):
                val = string.split(':')
                if val[0] not in acceptable_values:
                    ind2lose.append(index)
            junction_counts_df = junction_counts_df.append(junction_counts)
            junction_counts_df = junction_counts_df.drop(junction_counts.index[ind2lose])
            junction_counts_df = junction_counts_df[junction_counts_df['chromosome'].str.contains('chr')]
            junction_counts_df = junction_counts_df[junction_counts_df["chromosome"].str.contains("chrM") == False]
            junction_counts_df['intron_start'] = junction_counts_df['intron_start'] - 1
            junction_counts_df['intron_end'] = junction_counts_df['intron_end'] + 1
            junction_counts_df = junction_counts_df[junction_counts_df['intron_start'] != junction_counts_df['intron_end']]
            junc_string = []
            for index, row in junction_counts_df.iterrows():
                junc_string.append('{}:{}-{}'.format(row['chromosome'], row['intron_start'], row['intron_end']))
            junction_counts_df.insert(0, 'junction', junc_string)
            junction_counts_df = junction_counts_df.drop(['chromosome', 'intron_start', 'intron_end', 'strand', 'motif', 'is_annotated', 'multi_mapped_reads', 'max_overhang'], axis = 1)
            filename = file.split('SJ')[0]
            junction_counts_df.rename(columns = {'unique_reads':'{}'.format(filename)}, inplace = True)
            combined_junction = combined_junction.merge(junction_counts_df, on = 'junction', how = 'outer')
    
    combined_junction = combined_junction.fillna(0)
    combined_junction = combined_junction.reindex(sorted(combined_junction.columns), axis=1)
    junc_col = combined_junction['junction']
    combined_junction.drop(['junction'], axis=1, inplace = True)
    combined_junction.insert(0, 'junction', junc_col)
    return combined_junction

Z = combine_junction_counts(args.star_results_tumors)
Z.to_csv('{}/{}_junctions.txt'.format(args.output_directory, args.output_prefix), sep = '\t', index = False)

