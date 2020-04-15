# This script will generally follow suggestions in https://osca.bioconductor.org/index.html.
library(data.table)
library(SingleCellExperiment)

# Read it in and get in the correct format
if (!file.exists("data/celseq_matrix_ru10_molecules.tsv.725585")) {
  if (!file.exists("data/celseq_matrix_ru10_molecules.tsv.725585.gz")) {
    stop("can't find file!")
  }
  system("gunzip data/celseq_matrix_ru10_molecules.tsv.725585.gz")
}
umi_matrix <- data.table::fread("data/celseq_matrix_ru10_molecules.tsv.725585")
mat <- as.matrix(umi_matrix[,2:dim(umi_matrix)[2]])
rownames(mat) <- umi_matrix$gene
rm(umi_matrix)

sce <- SingleCellExperiment(assays = list(counts = mat))
# rm(mat)

# Get some metadata
map <- fread("data/celseq_meta.mapping2immport.725595.txt")
meta <- fread("data/SDY998-DR33_Subject_2_RNA_sequencing_result.txt")
meta_sub <- unique(meta[, .(`Subject Accession`, 
                        Species, 
                        Race, 
                        Gender, 
                        `Subject Age`, 
                        `ARM Accession`, 
                        `Biosample Accession`,
                        `Biosample Type`,
                        `Biosample Subtype`,
                        `Expsample Accession`)])
allmeta <- merge(map, meta_sub,
                 by.x = c("EXPSAMPLE_ACCESSION", "BIOSAMPLE_ACCESSION", "SUBJECT_ACCESSION"),
                 by.y = c("Expsample Accession", "Biosample Accession", "Subject Accession"))
rownames(allmeta) <- allmeta$cell_name

# What are all those cells missing from metadata???
# Just remove them for now
mat <- mat[,allmeta$cell_name]
sce <- SingleCellExperiment(assays = list(counts = mat),
                            colData = allmeta)
sce <- scater::addPerCellQC(sce)
sce[, sce$`ARM Accession` == "ARM3569"]

# Add some metadata about genes
sce <- scater::addPerFeatureQC(sce)

# Get gene ranges
# biomaRt
# annotationdbi
# Need a system to keep mapping updated -- first update gene name, then metadata about gene

# Get metadata from annotationhub
library(AnnotationHub)
ah <- AnnotationHub()
grch37 <- subset(ah, genome == "GRCh37")
grch37 <- grch37[["AH10684"]]

# Get ranges by gene name
ranges <- range(split(grch37, ~gene_name))

# TODO 
# There are a bunch of genes (134) which are missing from the ranges object... 
# For now I'll just remove them so I can move on. Come fix this later!
sce <- sce[rownames(sce) %in% names(ranges),]
rowRanges(sce) <- ranges[rownames(sce)]

