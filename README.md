# Single-cell transcriptomics identifies synergistic role of leukemic and non-leukemic immune repertoires in CD8+ T-cell large granular lymphocytic leukemia

Scripts to reproduce figures and analyses in the manuscript "Single-cell transcriptomics identifies synergistic role of leukemic and non-leukemic immune repertoires in CD8+ T-cell large granular lymphocytic leukemia" Huuhtanen, Bhattacharya et al., submitted 



# Pseudocode for CD8+ T-LGLL manuscript

## scTCRab-seq preprocessing

<pre>
<b>for</b> scTCRseq data, <i><b>do</b></i>  
  <b>read</b> 10X CellRanger output ## 11 T-LGLL and 6 healthy  
  <b>filter</b> ## based on the quality of cells (incomple TCRab information, loq confidence data)  
  <b>annotate</b> leukemic clones ## using manual annotation mark all clones as leukemic and non clones as possible non leukemic clones  
</pre>


## scRNAseq preprocessing

<pre>
<b>for</b> t_lgll_data, <i><b>do</b></i>  
  <b>read</b> 10X CellRanger output ## 11 T-LGLL samples  
  <b>filter</b> based on the quality of cells ## different for each individual to avoid loss of leukemic cells  

<b>for</b> healthy_data, <i><b>do</b></i>  
  <b>read</b> 10X CellRanger output 
  <b>filter</b> based on the quality of cells ## common for all healthy donors  

total_data = <b>merge</b> t_lgll_data_filtered with healthy_data_filtered
</pre>

<pre>
<b>for</b> total_data, <i><b>do</b></i> 
  <b>annotate</b> cells with singleR  
  <b>merge</b> with TCRdata  
  <b>calculate</b> latent represenation ## with scvi  
  <b>calculate</b> UMAP representation from latents  
  <b>calculate</b> clustering from latents  
  <b>decide</b> optimal clustering ## based on visual inspection  
  <b>scale</b> data  
</pre>

<pre>
<b>for</b> total_data, <i><b>do</b></i>  
  <b>calculate</b> DE-genes, pathways between clusters  
  <b>calculate</b> DE-genes, pathways between T-LGLL and healthy  
  <b>calculate</b> ligand-receptor pairs with cellphonedb  
  <b>calculate</b> regulons with scenic  
  <b>detect</b> STAT3 SNPs with vartrix  
</pre>
	
	
## Determine antigen-drive

<pre>
t_lgll_tcr_data = list scTCRab-seq, bulk-TCRb-seq samples

<b>for</b> t_lgll_tcr_data, <i><b>do</b></i>  
  <b>subsample</b> to 30k reads  
  <b>run</b> GLIPH2  
  <b>identify</b> t_lgll_clone  
  <i><b>if</b></i>   t_lgll_clone <i><b>in</b></i>   GLIPH2_results;  
    antigen_drive  
  <i><b>else</b></i>    
    no_antgen_drive  
</pre>

<pre>
other_tcr_data = list bulk-TCRb-seq from skcm, ra, healthy samples

<b>for</b> other_tcr_data, <i><b>do</b></i>
  <b>subsample</b> to 30k reads
  <b>run</b> GLIPH2
  <b>identify</b> largest_clone
  <i><b>if</b></i>   largest <i><b>in</b></i>   GLIPH2_results;  
    antigen_drive
  <i><b>else</b></i>    
    no_antgen_drive
</pre>

			 
	




## To reproduce the results:

### 1) Clone this repository

```
git clone https://github.com/janihuuh/cd8_tlgll_manu
cd path/to/cd8_tlgll_manu/
```

### 2) Obtain the data

* The processed scRNA+TCRab-seq data can be received from EGA (accession number EGAS00001005297, submission ongoing). 
* The processed bulk-RNA-seq data can be received from ArrayExpress (submission ongoing). 
* The TCRb-seq data can be received from immuneAccess (submission ongoing).

### 3) Create seurat-objects

```
Rscript R/main.R ## init the helper-functions, coloring, etc.
Rscript R/helper/run_preprocessTCRab.R ## preprocess the TCRab-seq
python python/run_scvi_example.py ## obtain the latent embeddings to use in clustering, UMAPs. This is just an example script
Rscript R/rnaseq/run_createSeurat.R ## read the CellRanger output, filter, merge TCRab-data, merge scVI latent embeddings, run singleR, DEs, pathways, etc.

```

### 4) Run additional RNAseq analyses

```
Rscript R/rnaseq/run_cellphone.R ## init data for CellPhoneDb
bash bash/run_cellphonedb.sh ## Run CellPhoneDb
Rscript R/rnaseq/run_scenic.R ## run Scenic-analysis
bash bash/run_vartrix.sh ## Run Vartrix
Rscript R/rnaseq/run_edgeR.R ## run bulk-RNAseq analysis
```

### 5) Run additional TCRseq analyses 

```
bash bash/run_vdjtools.sh ## preprocess the bulk-TCRb-data, subsample, calcualte diverisities, etc.
bash bash/run_gliph2_expample.sh ## run GLIPH2 collectively and on individual samples; this is just an example script
python python/run_tcrgp.py ## run TCRGP collectively on samples
Rscript R/tcrseq/run_antigen_drive.R ## run antigen-drive based analyses
```
