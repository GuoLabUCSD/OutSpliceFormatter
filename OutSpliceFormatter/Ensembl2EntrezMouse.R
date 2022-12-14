library(argparser)
suppressPackageStartupMessages(library(AnnotationDbi))
library(org.Mm.eg.db)
suppressPackageStartupMessages(library(dplyr))

parser = arg_parser('')

parser = add_argument(parser, '--output_directory', help = 'Output Directory Containing RSEM Expected Counts with ENSEMBL gene IDs')
parser = add_argument(parser, '--output_prefix', help = 'Filename Prefix')

args = parse_args(parser)

#Get IDs and Gene Symbols
#Remove Gene Copies on Pseudoautosomal Regions
my_rsem_data <- read.csv(sprintf('%s/%s_expected_counts_TMP.txt', args$output_directory, args$output_prefix), sep = '\t', header = TRUE)

my_rsem_data_clean <- my_rsem_data

my_rsem_data_clean <- my_rsem_data_clean[!grepl("PAR", my_rsem_data_clean$gene_id),]

my_rsem_data_clean$gene_id <- gsub('\\..*', '', my_rsem_data_clean$gene_id)

my_rsem_data_clean$ENTREZID <- mapIds(org.Mm.eg.db, keys=my_rsem_data_clean$gene_id, column="ENTREZID", keytype = "ENSEMBL", multiVals="first")

#my_rsem_data_clean$symbol <- mapIds(org.Mm.eg.db, keys=my_rsem_data_clean$gene_id, column="SYMBOL", keytype = "ENSEMBL", multiVals="first")

my_rsem_data_clean <- na.omit(my_rsem_data_clean)

#Match Ensembl IDs to corresponding Entrez IDs
my_rsem_data_clean <- my_rsem_data_clean[order(my_rsem_data_clean$gene_id),]
#my_rsem_data_clean <- my_rsem_data_clean %>% distinct(symbol, .keep_all = TRUE)
my_rsem_data_clean <- my_rsem_data_clean %>% distinct(ENTREZID, .keep_all = TRUE)

#my_rsem_data_clean$symbol_entrez <- paste(my_rsem_data_clean$symbol, my_rsem_data_clean$ENTREZID, sep = "|")

my_rsem_data_clean <- my_rsem_data_clean %>% dplyr::select(-gene_id)

my_rsem_data_final <- my_rsem_data_clean %>% dplyr::select(ENTREZID, everything())

colnames(my_rsem_data_final) <- gsub('\\.', "-", colnames(my_rsem_data_final))

#my_rsem_data_final <- my_rsem_data_final[order(my_rsem_data_final$symbol_entrez),]

#Save the Final Table
write.table(my_rsem_data_final, sprintf('%s/%s_expected_counts_entrez.txt', args$output_directory, args$output_prefix), sep = '\t', quote = FALSE, row.names = FALSE)

