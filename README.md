# exploreScRNASeq

I'm using this repo for some initial explorations into processing single-cell RNA seq data. 

Data comes from [ImmPort SDY998](https://www.immport.org/shared/study/SDY998).  
It can be downloaded through aspera from https://browser.immport.org/browser?path=SDY998

The scripts here expect the contents of `ResultFiles/RNA_sequencing_result` to be in a 
directory named `data`. Note that some of these files, when unzipped and loaded into R,
can take up mulitple gigabytes of memory. 


## SDY998

[https://www.immport.org/shared/study/SDY998](https://www.immport.org/shared/study/SDY998) 

- start date 2016
- 4 pubmed publications (2018-2019)
    - [Methods for high-dimensional analysis of cells dissociated from cryopreserved synovial tissue](https://www.ncbi.nlm.nih.gov/pubmed/29996944?dopt=Abstract)
    - [Defining inflammatory cell states in rheumatoid arthritis joint synovial tissues by integrating single-cell transcriptomics and mass cytometry.](https://www.ncbi.nlm.nih.gov/pubmed/31061532?dopt=Abstract)
    - [PD-1hiCXCR5- T peripheral helper cells promote B cell responses in lupus via MAF and IL-21.](https://www.ncbi.nlm.nih.gov/pubmed/31536480?dopt=Abstract)
    - [HBEGF+ macrophages identified in rheumatoid arthritis promote joint tissue invasiveness and are reshaped differentially by medications](https://www.biorxiv.org/content/10.1101/525758v1)

> The primary goal for RA arthroplasty P1 studies are: To establish if molecular signatures and pathways identified using core AMP technologies differ between OA and RA in 20 RA surgical samples and 10 OA arthroplasty samples.

AMP = [Accelerating medicines partnership](https://www.niams.nih.gov/grants-funding/funded-research/accelerating-medicines)

Data visualized here (via plotly) [https://singlecell.broadinstitute.org/single_cell/study/SCP279/amp-phase-1#study-summary](https://singlecell.broadinstitute.org/single_cell/study/SCP279/amp-phase-1#study-summary)

Available files

    data/
    ├── ReadMe_expression_2017-05-12.725584.docx
    ├── celseq_bad_barcodes.tsv.725588.gz
    ├── celseq_flow.tsv.725589.gz
    ├── celseq_flow_markers.tsv.725593.gz
    ├── celseq_matrix_ru10_molecules.tsv.725585.gz
    ├── celseq_matrix_ru10_reads.tsv.725590.gz
    ├── celseq_matrix_ru1_molecules.tsv.725583.gz
    ├── celseq_matrix_ru1_reads.tsv.725586.gz
    ├── celseq_meta.mapping2immport.725595.txt
    ├── celseq_meta.tsv.725591.gz
    ├── celseq_meta_unfiltered.tsv.725587.gz
    ├── celseq_molecule_counts.tsv.723951.gz
    ├── celseq_molecule_counts.tsv.725592.gz
    ├── celseq_star_log.tsv.725582.gz
    ├── low_input_gene_sample_tpm_matrix.725714.tsv
    ├── low_input_meta.725715.tsv

`celseq_molecule_counts` are the raw data — matrices generated from these

`_matrix_` files are GE matrices: each row is a gene, each column is a cell

`_molecules_` gives counts of UMI's

`_reads_` gives counts of reads 

`_ru10_` only includes molecules with ≥ 10 reads per UMI

`_ru1_` does not filter by ≥ 10 reads per UMI

`celseq_meta` metadata for all cells (see README for more info)
