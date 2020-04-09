# This script will generally follow suggestions in https://osca.bioconductor.org/index.html.


library(SingleCellExperiment)

# Read it in and get in the correct format
umi_matrix <- data.table::fread("data/celseq_matrix_ru10_molecules.tsv.725585.gz")
mat <- as.matrix(umi_matrix[,2:dim(umi_matrix)[2]])
rownames(mat) <- umi_matrix$gene
rm(umi_matrix)

sce <- SingleCellExperiment(assays = list(counts = mat))

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

sce <- SingleCellExperiment(assays = list(counts = mat),
                            colData = allmeta)

# What are all those cells missing from metadata???
# Just remove them for now
mat <- mat[,allmeta$cell_name]
sce <- SingleCellExperiment(assays = list(counts = mat),
                            colData = allmeta)
sce <- scater::addPerCellQC(sce)
sce[, sce$`ARM Accession` == "ARM3569"]

# Add some metadata about genes
sce <- scater::addPerFeatureQC(sce)
