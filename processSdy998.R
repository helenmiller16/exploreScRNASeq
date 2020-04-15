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
mat[is.na(mat)] <- 0
rownames(mat) <- umi_matrix$gene
rm(umi_matrix)

sce <- SingleCellExperiment(assays = list(counts = mat))
# rm(mat)

# Get some metadata
if (FALSE) {
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
  sce[, sce$`ARM Accession` == "ARM3569"]
}

# ----- QC -----
# Identify motochondrial transcripts 
# https://osca.bioconductor.org/quality-control.html#choice-of-qc-metrics

# NOTE: Need a system to keep mapping updated -- first update gene name, then metadata about gene


# Derived meta -----
# Get metadata from annotationhub
library(AnnotationHub)
# ah <- AnnotationHub()
# grch37 <- subset(ah, genome == "GRCh37")
# grch37 <- grch37[["AH10684"]]
grch37 <- AnnotationHub()[["AH10684"]]

# Get ranges by gene name
ranges <- range(split(grch37, ~gene_name))

# TODO 
# There are a bunch of genes (134) which are missing from the ranges object... 
# For now I'll just remove them so I can move on. Come fix this later!
sce <- sce[rownames(sce) %in% names(ranges),]
rowRanges(sce) <- ranges[rownames(sce)]

# Get percent mitochondrial dna
is.mito <- any(seqnames(rowRanges(sce))=="MT")

library(scater)
df <- perCellQCMetrics(sce, subsets=list(Mito=is.mito))
# sce <- addPerCellQC(sce, subsets=list(Mito=is.mito))

# Provided metadata ----- 
cell_meta <- fread("data/celseq_meta.tsv.725591.gz")
# rownames(cell_meta) <- cell_meta$cell_name
# colData(sce) <- cell_meta

# Remove low-quality cells 
# https://osca.bioconductor.org/quality-control.html#quality-control-outlier

# This is the same as quickPerCellQC(df, percent_subsets=c("subsets_Mito_percent"))
# This should be done separately on different batches


# QC by batch

cell_meta[, batch := paste0(sample, type)]
qc.mol <- isOutlier(cell_meta$molecules, log=TRUE, type="lower", batch = cell_meta$batch)
qc.nexprs <- isOutlier(cell_meta$genes_detected, log=TRUE, type="lower", batch = cell_meta$batch)
qc.mito <- isOutlier(cell_meta$percent_mt_molecules, type="higher", batch = cell_meta$batch)
discard <- qc.mol | qc.nexprs | qc.mito
DataFrame(LibSize=sum(qc.mol), NExprs=sum(qc.nexprs),
          MitoProp=sum(qc.mito), Total=sum(discard))


# Diagnostic plots 
# https://osca.bioconductor.org/quality-control.html#quality-control-plots

colData(sce) <- DataFrame(cell_meta)
sce$discard <- discard
plotColData(sce, x = "type", y = "molecules", colour_by = "discard", other_fields = "disease") + 
  facet_wrap(~disease) + 
  scale_y_log10() + ggtitle("Total count")
plotColData(sce, x = "type", y = "genes_detected", colour_by = "discard", other_fields = "disease") + 
  facet_wrap(~disease) + 
  scale_y_log10() + ggtitle("Detected Features")
plotColData(sce, x = "type", y = "percent_mt_molecules", colour_by = "discard", other_fields = "disease") + 
  facet_wrap(~disease) + 
  scale_y_log10() + ggtitle("Mito Percent")


plotColData(sce, x="molecules", y="percent_mt_molecules", 
            colour_by="discard", other_fields=c("disease", "type")) +
  facet_grid(type~disease) +
  theme(panel.border = element_rect(color = "grey"))

# Cell calling
# https://osca.bioconductor.org/quality-control.html#qc-droplets
if (FALSE) {
  library(DropletUtils)  
  library(Matrix)
  bcrank <- barcodeRanks(counts(sce))
  uniq <- !duplicated(bcrank$rank)
  plot(bcrank$rank[uniq], bcrank$total[uniq], log="xy",
       xlab="Rank", ylab="Total UMI count", cex.lab=1.2)
  
  set.seed(100)
  e.out <- emptyDrops(counts(sce))
  summary(e.out$FDR <= 0.001)
  table(Sig=e.out$FDR <= 0.001, Limited=e.out$Limited)
  limit <- 100   
  all.out <- emptyDrops(counts(sce), lower=limit, test.ambient=TRUE)
  hist(all.out$PValue[all.out$Total <= limit & all.out$Total > 0],
       xlab="P-value", main="", col="grey80") 
  
  sce <- sce.pbmc[,which(e.out$FDR <= 0.001)]
}


# Diagnose cell-type loss
# https://osca.bioconductor.org/quality-control.html#quality-control-discarded
discard <- sce$discard
filtered <- sce[, !discard]
lost <- calculateAverage(counts(sce)[,!discard])
kept <- calculateAverage(counts(sce)[,discard])
library(edgeR)
logged <- cpm(cbind(lost, kept), log=TRUE, prior.count=2)
logFC <- logged[,1] - logged[,2]
abundance <- rowMeans(logged)
plot(abundance, logFC, xlab="Average count", ylab="Log-FC (lost/kept)", pch=16)
points(abundance[is.mito], logFC[is.mito], col="dodgerblue", pch=16)
points(abundance[is.mito], logFC[is.mito], col="dodgerblue", pch=16)
