# Read raw counts
library(data.table)

# WARNING!!!!! You need ~ 9 gb of memory to read this in!!!!
# It takes up < 7 gb once loaded (>150m rows)
if (FALSE) {
  raw <- fread("data/celseq_molecule_counts.tsv.723951.gz")
}
# raw1 <- fread("data/celseq_molecule_counts.tsv.725592.gz")

# ------ Unfiltered read counts ------
# This guy still takes up >1gb!
reads_matrix <- fread("data/celseq_matrix_ru1_reads.tsv.725586.gz")
dim(reads_matrix)
reads_matrix[1:10, 1:10]
# Lots of NA! Presumably that means 0. What portion of just the first 10 cells are NOT missing? 
sum(!is.na(reads_matrix[,1:10]))/(10*dim(reads_matrix)[1])

# Let's look at some more things from the first 10 cells
sub <- reads_matrix[,1:10]
range(unlist(sub[,2:10]), na.rm = TRUE)
hist(unlist(sub[,2:10]), 100)
hist(log(unlist(sub[,2:10])), 100)

# Get the mean of each cell (should be about the same)
means <- apply(sub[,2:10], 2, mean, na.rm = TRUE)
# That actually seems like a pretty big range -- but maybe I need to set NA's to 0. 

rm(reads_matrix)
rm(sub)
# ------ Filtered UMI counts ------
umi_matrix <- fread("data/celseq_matrix_ru10_molecules.tsv.725585.gz")
dim(umi_matrix)
umi_matrix[1:10, 1:10]
# a bunch of genes are filtered out

sum(!is.na(umi_matrix[,1:10]))/(10*dim(umi_matrix)[1])
# 13.79% not missing

sub <- umi_matrix[,1:10]
range(unlist(sub[,2:10]), na.rm = TRUE)
# Range is WAY smaller than total reads!
hist(unlist(sub[,2:10]), 100)
hist(log(unlist(sub[,2:10])), 100)
means <- apply(sub[,2:10], 2, mean, na.rm = TRUE)
means
hist(means)
# still seems like a pretty big range...

# ------ low_input_gene_sample_tpm_matrix -----
# Not sure what this is...
tpm <- fread("data/low_input_gene_sample_tpm_matrix.725714.tsv")
dim(tpm)
names(tpm)
meta <- fread("data/low_input_meta.725715.tsv")
meta


# ------ mapping2immport ------
map <- fread("data/celseq_meta.mapping2immport.725595.txt")
map

# ----- more metadata ------
meta <- fread("data/SDY998-DR33_Subject_2_RNA_sequencing_result.txt")
meta
# Lost of good metadata in here!
# this looks like what might go into database...? 
