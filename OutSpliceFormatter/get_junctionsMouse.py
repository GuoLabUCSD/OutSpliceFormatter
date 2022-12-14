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

def fix_mouse_chrom(junction_counts):
    chrom_dict = {'1': 'chr1', '2': 'chr2', '3': 'chr3', '4': 'chr4',
                  '5': 'chr5', '6': 'chr6', '7': 'chr7', '8': 'chr8',
                  '9': 'chr9', '10': 'chr10', '11': 'chr11', '12': 'chr12',
                  '13': 'chr13', '14': 'chr14', '15': 'chr15', '16': 'chr16',
                  '17': 'chr17', '18': 'chr18', '19': 'chr19', 'X': 'chrX', 
                  'Y': 'chrY'}
    junction_counts = junction_counts.replace({'chromosome': chrom_dict})
    return junction_counts

def combine_junction_counts(star_results_tumors, star_results_normals):
    combined_junction = pd.DataFrame(columns = ['junction'])

    for file in os.listdir(star_results_tumors):
        junction_counts_df = pd.DataFrame()
        if file.endswith("SJ.out.tab"):
            path = os.path.join(star_results_tumors, file)
            junction_counts = pd.read_csv(path, sep = '\t', header = None, names = ['chromosome', 'intron_start', 'intron_end', 'strand', 'motif', 'is_annotated', 'unique_reads', 'multi_mapped_reads', 'max_overhang'], dtype={'chromosome': 'str'})
            junction_counts = fix_mouse_chrom(junction_counts)
            junction_counts_df = junction_counts_df.append(junction_counts)
            junction_counts_df = junction_counts_df[junction_counts_df['chromosome'].str.contains('chr')]
            junction_counts_df = junction_counts_df[junction_counts_df["chromosome"].str.contains("chrM") == False]
            junction_counts_df['intron_start'] = junction_counts_df['intron_start'] - 1
            junction_counts_df['intron_end'] = junction_counts_df['intron_end'] + 1
            junction_counts_df = junction_counts_df[junction_counts_df['intron_start'] != junction_counts_df['intron_end']]
            junc_string = []
            for index, row in junction_counts_df.iterrows():
                junc_string.append('{}:{}-{}'.format(row['chromosome'], row['intron_start'], row['intron_end']))
                #if row['strand'] == 1:
                #    junc_string.append('{}:{}:+,{}:{}:+'.format(row['chromosome'], row['intron_start'], row['chromosome'], row['intron_end']))
                #elif row['strand'] == 2:
                #    junc_string.append('{}:{}:-,{}:{}:-'.format(row['chromosome'], row['intron_start'], row['chromosome'], row['intron_end']))
                #elif row['strand'] == 0:
                #    junc_string.append('{}:{}:?,{}:{}:?'.format(row['chromosome'], row['intron_start'], row['chromosome'], row['intron_end']))
            junction_counts_df.insert(0, 'junction', junc_string)
            junction_counts_df = junction_counts_df.drop(['chromosome', 'intron_start', 'intron_end', 'strand', 'motif', 'is_annotated', 'multi_mapped_reads', 'max_overhang'], axis = 1)
            filename = file.split('SJ')[0]
            #filename = filename.ljust(10, 'X')
            junction_counts_df.rename(columns = {'unique_reads':'{}'.format(filename)}, inplace = True)
            combined_junction = combined_junction.merge(junction_counts_df, on = 'junction', how = 'outer')
            
    for file in os.listdir(star_results_normals):
        junction_counts_df = pd.DataFrame()
        if file.endswith("SJ.out.tab"):
            path = os.path.join(star_results_normals, file)
            junction_counts = pd.read_csv(path, sep = '\t', header = None, names = ['chromosome', 'intron_start', 'intron_end', 'strand', 'motif', 'is_annotated', 'unique_reads', 'multi_mapped_reads', 'max_overhang'], dtype={'chromosome': 'str'})
            junction_counts = fix_mouse_chrom(junction_counts)
            junction_counts_df = junction_counts_df.append(junction_counts)
            junction_counts_df = junction_counts_df[junction_counts_df['chromosome'].str.contains('chr')]
            junction_counts_df = junction_counts_df[junction_counts_df["chromosome"].str.contains("chrM") == False]
            junction_counts_df['intron_start'] = junction_counts_df['intron_start'] - 1
            junction_counts_df['intron_end'] = junction_counts_df['intron_end'] + 1
            junction_counts_df = junction_counts_df[junction_counts_df['intron_start'] != junction_counts_df['intron_end']]
            junc_string = []
            for index, row in junction_counts_df.iterrows():
                junc_string.append('{}:{}-{}'.format(row['chromosome'], row['intron_start'], row['intron_end']))
                #if row['strand'] == 1:
                #    junc_string.append('{}:{}:+,{}:{}:+'.format(row['chromosome'], row['intron_start'], row['chromosome'], row['intron_end']))
                #elif row['strand'] == 2:
                #    junc_string.append('{}:{}:-,{}:{}:-'.format(row['chromosome'], row['intron_start'], row['chromosome'], row['intron_end']))
                #elif row['strand'] == 0:
                #    junc_string.append('{}:{}:?,{}:{}:?'.format(row['chromosome'], row['intron_start'], row['chromosome'], row['intron_end']))
            junction_counts_df.insert(0, 'junction', junc_string)
            junction_counts_df = junction_counts_df.drop(['chromosome', 'intron_start', 'intron_end', 'strand', 'motif', 'is_annotated', 'multi_mapped_reads', 'max_overhang'], axis = 1)
            filename = file.split('SJ')[0]
            #filename = filename.ljust(10, 'X')
            junction_counts_df.rename(columns = {'unique_reads':'{}'.format(filename)}, inplace = True)
            combined_junction = combined_junction.merge(junction_counts_df, on = 'junction', how = 'outer')
    
    combined_junction = combined_junction.fillna(0)
    combined_junction = combined_junction.reindex(sorted(combined_junction.columns), axis=1)
    junc_col = combined_junction['junction']
    combined_junction.drop(['junction'], axis=1, inplace = True)
    combined_junction.insert(0, 'junction', junc_col)
    #combined_junction.loc[-1] = 'raw_counts'
    #combined_junction.index = combined_junction.index + 1
    #combined_junction = combined_junction.sort_index()
    return combined_junction

Z = combine_junction_counts(args.star_results_tumors, args.star_results_normals)
Z.to_csv('{}/{}_junctions.txt'.format(args.output_directory, args.output_prefix), sep = '\t', index = False)
