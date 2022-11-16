# OutSpliceFormatter:
This tool aims to ease a user's RNA-Seq Splicing Analyses with the OutSplice Software on Human or Mouse Genomes by automatically extracting the total and junction read counts from a STAR alignment output directory. It will then output a matrix with all of the data combined together that can be easily provided to OutSplice for immediate use. Optionally, if RSEM directories are provided, the tool will also extract the expected counts expression results and perform a quartile normalization based on the TCGA mRNA-seq pipeline, outputting a normalized expression matrix for use with OutSplice. If the user wishes to use their own normalization method, a matrix of all expected counts is also provided.

## Pre-Requisites:
 * STAR Results for Tumor and Normal Samples
	* For human data: Please ensure the chromsome column in your SJ.out.tab file has the below labelling convention
				Ex) chr1, chr2, etc.
	* For mouse data: Please ensure the chromsome column in your SJ.out.tab file has the below labelling convention
				Ex) 1, 2, etc.
 * RSEM Results for Tumor and Normal Samples [OPTIONAL]
	

## Installation:
	1. git clone (https://github.com/GuoLabUCSD/OutSplice)
	2. Download quartile_norm.pl (https://github.com/mozack/ubu/blob/master/src/perl/quartile_norm.pl)
	3. Place the quartile_norm.pl file inside the OutSplice_Formatter directory


## Minimum Requirements:
 * Python 3 with the following modules installed:
	* pandas (>= version 1.4.3)


## Requirements for RSEM Data Formatting Only:
 * Perl (>= version 5.26.3)
	* quartile_norm.pl
 * R (>= version 4.1)
	* AnnotationDbi (>= version 1.56.2)
	* argparser (>= version 0.7.1)
	* org.Hs.eg.db AND/OR org.Mm.eg.db (>= version 3.14.0)


## Arguments:

	-a      
		Path to a Directory containing STAR's Log.final.out and SJ.tab.out files for the tumor samples
	-b      
		Path to a Directory containing STAR's Log.final.out and SJ.tab.out files for the normal samples
	-c      
		Path to a Directory containing RSEM's genes.results files for the tumor samples [OPTIONAL]
	-d      
		Path to a Directory containing RSEM's genes.results files for the normal samples [OPTIONAL]
	-o      
		Output Directory to store files
	-p      
		Prefix to what the output files should be named
	-m      
		Data is from a Mouse Genome [OPTIONAL]

## Example Usage:
* From the OutSpliceFormatter Directory

		./format_data.sh -a ../STAR_OUTPUT_TUMORS/ -b ../STAR_OUTPUT_NORMALS -c ../RSEM_OUTPUT_TUMORS -d ../RSEM_OUTPUT_NORMALS -o matrix_files/ -p my_samples -m

## Output:
	{Filename Prefix}_phenos.txt:
		Tab separated text file containing a matrix of sample names and whether they are tumors (True) or normals (False)
	
	{Filename Prefix}_rawcounts.txt:
		Tab separated text file containing a matrix of the number of uniquely mapped reads for each sample based on STAR's Log.final.out file

	{Filename Prefix}_expected_counts_entrez.txt:
		Tab separated text file containing a matrix of the expected counts expression results from RSEM's genes.results file with entrez IDs

	{Filename Prefix}_genes_normalized.txt:
		Tab separated text file containing a matrix of the quartile normalized expression results from RSEM's genes.results file with entrez IDs

	{Filename Prefix}_junctions.txt:
		Tab separated text file containing a matrix of junction read counts based on STAR's SJ.out.tab file
