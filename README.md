# Single-cell transcriptomics identifies synergistic role of leukemic and non-leukemic immune repertoires in CD8+ T-cell large granular lymphocytic leukemia

Scripts to reproduce figures and analyses in the manuscript "Single-cell transcriptomics identifies synergistic role of leukemic and non-leukemic immune repertoires in CD8+ T-cell large granular lymphocytic leukemia" Huuhtanen, Bhattacharya et al., submitted 

## To reproduce the results:

### First, clone this repository

```
git clone https://github.com/janihuuh/cd8_tlgll_manu
cd path/to/cd8_tlgll_manu/
```

#### Obtain the data

* The processed scRNA+TCRab-seq data can be received from EGA (accession number EGAS00001005297, submission ongoing). 
* The processed bulk-RNA-seq data can be received from ArrayExpress (submission ongoing). 
* The TCRb-seq data can be received from immuneAccess (submission ongoing).

### To create seurat-objects, you need to run

```
Rscript R/main.R ## init the helper-functions, coloring, etc.
Rscript R/helper/run_preprocessTCRab.R ## preprocess the TCRab-seq
python python/run_scvi_example.py ## obtain the latent embeddings to use in clustering, UMAPs. This is just an example script
Rscript R/rnaseq/run_createSeurat.R ## read the CellRanger output, filter, merge TCRab-data, merge scVI latent embeddings, run singleR, DEs, pathways, etc.

```

### To analyze other scRNAseq aspects, you need to run

```
Rscript R/rnaseq/run_cellphone.R ## init data for CellPhoneDb
bash bash/run_cellphonedb.sh ## Run CellPhoneDb
Rscript R/rnaseq/run_scenic.R ## run Scenic-analysis
```

### To analyze TCR-seq data, you need to run

```
bash bash/run_vdjtools.sh ## preprocess the bulk-TCRb-data, subsample, calcualte diverisities, etc.
bash bash/run_gliph2_expample.sh ## run GLIPH2 collectively and on individual samples; this is just an example script
python python/run_tcrgp.py ## run TCRGP collectively on samples
Rscript R/tcrseq/run_antigen_drive.R ## run antigen-drive based analyses
```
