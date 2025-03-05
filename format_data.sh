#! /bin/bash

#Get required arguments
usage() {
  cat <<EOF
Usage: $(basename "${BASH_SOURCE[0]}") <args>

This script will extract output from your STAR alignment and RSEM Expression files and put them into the proper format for OutSplice.

Arguments:

-a 	Path to a Directory containing STAR's Log.final.out and SJ.tab.out files for the tumor samples
-b 	Path to a Directory containing STAR's Log.final.out and SJ.tab.out files for the normal samples
-c	Path to a Directory containing RSEM's genes.results files for the tumor samples [OPTIONAL]
-d	Path to a Directory containing RSEM's genes.results files for the normal samples [OPTIONAL]
-o 	Output Directory to store files
-p 	Prefix to what the output files should be named
-l	Path to the OutSplice Formatter Directory
-m	Data is from a Mouse Genome [OPTIONAL]
EOF
  exit
}

while getopts "a:b:c:d:o:p:l:mh" flag ; do

	case "${flag}" in
		a) star_results_tumors=${OPTARG};;
		b) star_results_normals=${OPTARG};;
		c) rsem_results_tumors=${OPTARG};;
		d) rsem_results_normals=${OPTARG};;
		o) output_directory=${OPTARG};;
		p) output_prefix=${OPTARG};;
		l) formatter_directory=${OPTARG};;
		m) mouse='This_is_a_Mouse_Genome';;
		h) usage;;
	esac
done

if [ -z $star_results_tumors ] || [ -z $star_results_normals ] || [ -z $output_directory ] || [ -z $output_prefix ] || [ -z $formatter_directory ]; then
	echo "Missing required arguments"
	usage
	exit 1
fi

if [ ! -d $output_directory ]; then
	echo "The specified output directory does not exist"
	exit 1
fi

#Get phenotype matrix
echo "Generating Phenotype Matrix"

python $formatter_directory/get_phenotypes.py \
--star_results_tumors $star_results_tumors \
--star_results_normals $star_results_normals \
--output_directory $output_directory \
--output_prefix $output_prefix

#Get rawcounts file
echo "Getting Rawcounts"

python $formatter_directory/get_rawcounts.py \
--star_results_tumors $star_results_tumors \
--star_results_normals $star_results_normals \
--output_directory $output_directory \
--output_prefix $output_prefix

echo "Rawcounts File Generated"

#Get Expected Counts from RSEM
if [ -z $rsem_results_tumors ] && [ -z $rsem_results_normals ]; then
	echo "No RSEM data provided, skipping data extraction..."
else
	echo "Getting Expecting Counts"

	python $formatter_directory/get_expectedCounts_RSEM.py --rsem_results_tumors $rsem_results_tumors --rsem_results_normals $rsem_results_normals --output_directory $output_directory --output_prefix $output_prefix

	echo "Converting IDs"
fi

#Convert ENSEMBL to ENTREZ
if [ -z $rsem_results_tumors ] && [ -z $rsem_results_normals ]; then
	echo "No RSEM data provided, skipping ID conversion..."
else
	if [ -z $mouse ]; then
		Rscript $formatter_directory/Ensembl2Entrez.R --output_directory $output_directory --output_prefix $output_prefix
	else
		Rscript $formatter_directory/Ensembl2EntrezMouse.R --output_directory $output_directory --output_prefix $output_prefix
	fi

	echo "IDs Converted"
fi


#RUN RSEM quartile normalization from the TCGA's Pipeline (Source: https://github.com/mozack/ubu/blob/master/src/perl/quartile_norm.pl)
if [ -z $rsem_results_tumors ] && [ -z $rsem_results_normals ]; then
	echo "Skipping Normalization"
else
	echo "Executing Quartile Normalization"

	perl $formatter_directory/quartile_norm.pl -s 1 -c -1 -q 75 -t 1000 -o $output_directory/$output_prefix\_genes_normalized_TMP.txt $output_directory/$output_prefix\_expected_counts_entrez.txt
fi


#Format Normalized Data
if [ -z $rsem_results_tumors ] && [ -z $rsem_results_normals ]; then
	echo "Skipping RSEM Formatting"
else
	echo 'Formatting'

	python $formatter_directory/format_normalized_RSEM.py --output_directory $output_directory --output_prefix $output_prefix

	echo 'Normalized Expression File Generated'
fi

#Get junction counts from STAR
echo 'Getting Junction Counts'

python $formatter_directory/get_junctions.py --star_results_tumors $star_results_tumors --star_results_normals $star_results_normals --output_directory $output_directory --output_prefix $output_prefix

echo "Junction's File Generated"

#Remove Temporary Files
if [ -z $rsem_results_tumors ] && [ -z $rsem_results_normals ]; then
	echo "All done!"
else
	rm $output_directory/$output_prefix\_expected_counts_TMP.txt
	rm $output_directory/$output_prefix\_genes_normalized_TMP.txt

	echo "All done!"
fi
 
