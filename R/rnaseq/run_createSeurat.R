
## Run create seurat files

## Fig 1b
lgll_dirs      <- list.dirs("data/scRNAseq/lgll/", recursive = T) %>% grep(pattern = "filtered_feature_bc_matrix", value = T)
healthy_dirs   <- list.dirs("data/scRNAseq/healthy/", recursive = T) %>% grep(pattern = "filtered_feature_bc_matrix", value = T)
folders        <- c(lgll_dirs, healthy_dirs)
scrnaseq_files <- lapply(folders, function(x){message(x); Read10X(data.dir = x) %>% CreateSeuratObject(project = extractSeuratName(x), min.cells = 3, min.features = 200)})
lgll_healthy_seurat  <- merge(scrnaseq_files[[1]], scrnaseq_files[-1], add.cell.ids = extractSeuratName(folders))

## Basic QC
lgll_healthy_seurat  <- PercentageFeatureSet(lgll_healthy_seurat, pattern = "^MT-", col.name = "percent.mt")
lgll_healthy_seurat  <- PercentageFeatureSet(lgll_healthy_seurat, pattern = "^RP", col.name = "percent.ribo")
lgll_healthy_seurat  <- PercentageFeatureSet(lgll_healthy_seurat, features = cycle.genes, col.name = "percent.cycle")
lgll_healthy_seurat@meta.data$barcode   <- colnames(lgll_healthy_seurat)
lgll_healthy_seurat  <- lgll_healthy_seurat %>% getQC()

## Put TCR
tot_barcode <- fread("data/scRNAseq+TCRseq/preprocessed/tot_healthy_barcode.txt")
lgll_healthy_seurat <- mergeTCRtoSeurat(seurat_object = lgll_healthy_seurat, tcr_df = tot_barcode)

## Get latents, use them for UMAPs, clustering
lgll_healthy_seurat %>% getScviInput()
latents <- fread("results/scvi/latents_lgll_healthy.txt")

## Now run python/run_scvi_example.py to obtain the latent dimensions
lgll_healthy_seurat <- lgll_healthy_seurat %>% putLatentsSeurat()
lgll_healthy_seurat <- lgll_healthy_seurat %>% getLatentUMAP() %>% getLatentClustering()

## Scale data
clonality_genes <- getClonalityGenes(lgll_healthy_seurat)
unwanted_genes  <- getUnwantedGenes(lgll_healthy_seurat)
lgll_healthy_seurat <- lgll_healthy_seurat %>% preprocessSeurat(cells.to.use = colnames(lgll_healthy_seurat))

Idents(lgll_healthy_seurat) <- lgll_healthy_seurat$RNA_snn_res.0.2
lgll_healthy_seurat$cluster <- Idents(lgll_healthy_seurat)

## Get singleR
lgll_healthy_seurat <- lgll_healthy_seurat %>% getSingler()

## Get DE-genes
lgll_healthy_cluster_markers <- FindAllMarkers(lgll_healthy_seurat, test.use = "t", max.cells.per.ident = 3e3) %>% filter(p_val_adj < 0.05 & ave_logFC > 0)
lgll_healthy_deg_markers     <- lapply(unique(lgll_healthy_seurat$cluster), getDEGbyClusterLGL)

lgl_up_hallmark <- lgl_healthy_markers %>% filter(avg_logFC > 0) %>% getHypergeometric(universe_df = universe_df, term_df = hallmark) %>% filter(p.adjust < 0.05)
lgl_dn_hallmark <- lgl_healthy_markers %>% filter(avg_logFC < 0) %>% getHypergeometric(universe_df = universe_df, term_df = hallmark) %>% filter(p.adjust < 0.05)
lgl_up_go       <- lgl_healthy_markers %>% filter(avg_logFC > 0) %>% getHypergeometric(universe_df = universe_df, term_df = go) %>% filter(p.adjust < 0.05)
lgl_dn_go       <- lgl_healthy_markers %>% filter(avg_logFC < 0) %>% getHypergeometric(universe_df = universe_df, term_df = go) %>% filter(p.adjust < 0.05)




## Fig 1c
cells.to.keep           <- lgll_healthy_seurat@meta.data %>% filter(!is.na(new_clontoypes_id)) %>% pull(barcode)
lgll_healthy_tcr_seurat <- subset(lgll_healthy_seurat, cells = cells.to.keep)
lgll_healthy_tcr_seurat <- lgll_healthy_tcr_seurat %>% getLatentUMAP() %>% getLatentClustering()
lgll_healthy_tcr_seurat <- lgll_healthy_tcr_seurat %>% preprocessSeurat(cells.to.use = colnames(lgll_healthy_seurat))

Idents(lgll_healthy_tcr_seurat) <- lgll_healthy_seurat$RNA_snn_res.0.5
lgll_healthy_tcr_seurat$cluster <- Idents(lgll_healthy_tcr_seurat)


lgll_healthy_tcr_cluster_markers <- FindAllMarkers(lgll_healthy_tcr_seurat, test.use = "t", max.cells.per.ident = 3e3) %>% filter(p_val_adj < 0.05 & ave_logFC > 0)
lgll_healthy_tcr_deg_markers     <- lapply(unique(lgll_healthy_tcr_seurat$cluster), getDEGbyClusterLGL)

lgl_tcr_up_hallmark <- lgll_healthy_tcr_deg_markers %>% filter(avg_logFC > 0) %>% getHypergeometric(universe_df = universe_df, term_df = hallmark) %>% filter(p.adjust < 0.05)
lgl_tcr_dn_hallmark <- lgll_healthy_tcr_deg_markers %>% filter(avg_logFC < 0) %>% getHypergeometric(universe_df = universe_df, term_df = hallmark) %>% filter(p.adjust < 0.05)
lgl_tcr_up_go       <- lgll_healthy_tcr_deg_markers %>% filter(avg_logFC > 0) %>% getHypergeometric(universe_df = universe_df, term_df = go) %>% filter(p.adjust < 0.05)
lgl_tcr_dn_go       <- lgll_healthy_tcr_deg_markers %>% filter(avg_logFC < 0) %>% getHypergeometric(universe_df = universe_df, term_df = go) %>% filter(p.adjust < 0.05)



## Fig 1d
clonotypes.to.keep               <- lgll_healthy_tcr_seurat@meta.data %>% group_by(new_clonotypes_id) %>% summarise(n=n()) %>% filter(n>=10) %>% pull(new_clonotypes_id)
cells.to.keep                    <- lgll_healthy_tcr_seurat@meta.data %>% filter(new_clontoypes_id %in% clonotypes.to.keep) %>% pull(barcode)
lgll_healthy_hyperexp_tcr_seurat <- subset(lgll_healthy_tcr_seurat, cells = cells.to.keep)
lgll_healthy_hyperexp_tcr_seurat <- lgll_healthy_tcr_seurat %>% lgll_healthy_hyperexp_tcr_seurat()
lgll_healthy_hyperexp_tcr_seurat <- lgll_healthy_hyperexp_tcr_seurat %>% preprocessSeurat(cells.to.use = colnames(lgll_healthy_seurat))

lgll_healthy_hyperexp_cluster_markers <- FindAllMarkers(lgll_healthy_hyperexp_seurat, test.use = "t", max.cells.per.ident = 3e3) %>% filter(p_val_adj < 0.05 & ave_logFC > 0)
lgll_healthy_hyperexp_deg_markers     <- lapply(unique(lgll_healthy_hyperexp_seurat$cluster), getDEGbyClusterLGL)

lgl_hyperexp_up_hallmark <- lgll_healthy_hyperexp_deg_markers %>% filter(avg_logFC > 0) %>% getHypergeometric(universe_df = universe_df, term_df = hallmark) %>% filter(p.adjust < 0.05)
lgl_hyperexp_dn_hallmark <- lgll_healthy_hyperexp_deg_markers %>% filter(avg_logFC < 0) %>% getHypergeometric(universe_df = universe_df, term_df = hallmark) %>% filter(p.adjust < 0.05)
lgl_hyperexp_up_go       <- lgll_healthy_hyperexp_deg_markers %>% filter(avg_logFC > 0) %>% getHypergeometric(universe_df = universe_df, term_df = go) %>% filter(p.adjust < 0.05)
lgl_hyperexp_dn_go       <- lgll_healthy_hyperexp_deg_markers %>% filter(avg_logFC < 0) %>% getHypergeometric(universe_df = universe_df, term_df = go) %>% filter(p.adjust < 0.05)





## Fig 2b
lgll_dirs      <- list.dirs("data/scRNAseq/lgll/", recursive = T) %>% grep(pattern = "filtered_feature_bc_matrix", value = T)
scrnaseq_files <- lapply(lgll_dirs, function(x){message(x); Read10X(data.dir = x) %>% CreateSeuratObject(project = extractSeuratName(x), min.cells = 3, min.features = 200)})
lgll_seurat    <- merge(scrnaseq_files[[1]], scrnaseq_files[-1], add.cell.ids = extractSeuratName(folders))

## Basic QC
lgll_seurat  <- PercentageFeatureSet(lgll_seurat, pattern = "^MT-", col.name = "percent.mt")
lgll_seurat  <- PercentageFeatureSet(lgll_seurat, pattern = "^RP", col.name = "percent.ribo")
lgll_seurat  <- PercentageFeatureSet(lgll_seurat, features = cycle.genes, col.name = "percent.cycle")
lgll_seurat@meta.data$barcode <- colnames(lgll_seurat)
lgll_seurat  <- lgll_seurat %>% getLgllQC()
lgll_seurat  <- mergeTCRtoSeurat(seurat_object = lgll_seurat, tcr_df = tot_barcode)


## Now run python/run_scvi_example.py to obtain the latent dimensions
lgll_clones_seurat %>% getScviInput()
latents <- fread("results/scvi/latents_lgll_clones.txt")

## Get latents, use them for UMAPs, clustering
lgll_clones_seurat <- lgll_clones_seurat %>% putLatentsSeurat()
lgll_clones_seurat <- lgll_clones_seurat %>% getLatentUMAP() %>% getLatentClustering()
lgll_clones_seurat <- lgll_clones_seurat %>% preprocessSeurat(cells.to.use = colnames(lgll_clones_seurat))

Idents(lgll_clones_seurat) <- lgll_clones_seurat$RNA_snn_res.0.2
lgll_clones_seurat$cluster <- Idents(lgll_clones_seurat)

lgll_clones_markers <- FindAllMarkers(lgll_clones_seurat, test.use = "t", max.cells.per.ident = 3e3) %>% filter(p_val_adj < 0.05 & ave_logFC > 0)







## Fig 4a
lgll_dirs    <- list.dirs("data/scRNAseq/lgll/", recursive = T) %>% grep(pattern = "filtered_feature_bc_matrix", value = T)
healthy_dirs <- list.dirs("data/scRNAseq/healthy/", recursive = T) %>% grep(pattern = "filtered_feature_bc_matrix", value = T)
cml_dirs     <- list.dirs("data/scRNAseq/cml/", recursive = T) %>% grep(pattern = "filtered_feature_bc_matrix", value = T)
cll_dirs     <- list.dirs("data/scRNAseq/cll/", recursive = T) %>% grep(pattern = "filtered_feature_bc_matrix", value = T)
nsclc_dirs   <- list.dirs("data/scRNAseq/nsclc/", recursive = T) %>% grep(pattern = "filtered_feature_bc_matrix", value = T)
rcc_dirs     <- list.dirs("data/scRNAseq/rcc/", recursive = T) %>% grep(pattern = "filtered_feature_bc_matrix", value = T)

folders              <- c(lgll_dirs, healthy_dirs, cml_dirs, cll_dirs, nsclc_dirs, rcc_dirs)
scrnaseq_files       <- lapply(folders, function(x){message(x); Read10X(data.dir = x) %>% CreateSeuratObject(project = extractSeuratName(x), min.cells = 3, min.features = 200)})
non_leukemic_seurat  <- merge(scrnaseq_files[[1]], scrnaseq_files[-1], add.cell.ids = extractSeuratName(folders))

## Basic QC
non_leukemic_seurat  <- PercentageFeatureSet(non_leukemic_seurat, pattern = "^MT-", col.name = "percent.mt")
non_leukemic_seurat  <- PercentageFeatureSet(non_leukemic_seurat, pattern = "^RP", col.name = "percent.ribo")
non_leukemic_seurat  <- PercentageFeatureSet(non_leukemic_seurat, features = cycle.genes, col.name = "percent.cycle")
non_leukemic_seurat@meta.data$barcode   <- colnames(non_leukemic_seurat)
non_leukemic_seurat  <- non_leukemic_seurat %>% getQC()

## Put TCR
tot_barcode <- fread("data/scRNAseq+TCRseq/preprocessed/tot_healthy_barcode.txt")
non_leukemic_seurat <- mergeTCRtoSeurat(seurat_object = non_leukemic_seurat, tcr_df = tot_barcode)

## Remove leukemic cells
lgll.cells <- colnames(lgll_seurat)
non_leukemic_seurat <- subset(non_leukemic_seurat, lgll.cells)

## Now run python/run_scvi_example.py to obtain the latent dimensions
non_leukemic_seurat %>% getScviInput()
latents <- fread("results/scvi/latents_lgll_healthy.txt")

## Get latents, use them for UMAPs, clustering
non_leukemic_seurat <- non_leukemic_seurat %>% putLatentsSeurat()
non_leukemic_seurat <- non_leukemic_seurat %>% getLatentUMAP() %>% getLatentClustering()

## Get rid of CLL, CML and healthy-specific (batch) clusters
cells.to.keep       <- non_leukemic_seurat@meta.data %>% filter(RNA_snn_res.0.5 %in% c(5,6,27)) %>% pull(barcode)
non_leukemic_seurat <- subset(non_leukemic_seurat, cells.to.keep)
non_leukemic_seurat <- non_leukemic_seurat %>% getLatentUMAP() %>% getLatentClustering()

## Scale data
clonality_genes <- getClonalityGenes(non_leukemic_seurat)
unwanted_genes  <- getUnwantedGenes(non_leukemic_seurat)
non_leukemic_seurat <- non_leukemic_seurat %>% preprocessSeurat(cells.to.use = colnames(non_leukemic_seurat))

## Get singleR
non_leukemic_seurat <- non_leukemic_seurat %>% getSingler()

Idents(non_leukemic_seurat) <- non_leukemic_seurat$RNA_snn_res.0.5
non_leukemic_seurat$cluster <- Idents(non_leukemic_seurat)

## Get DE-genes and pathways
non_leukemic_cluster_markers <- FindAllMarkers(non_leukemic_seurat, test.use = "t", max.cells.per.ident = 3e3) %>% filter(p_val_adj < 0.05 & ave_logFC > 0)


## a) T-LGLL to healthy
cells.to.keep               <- non_leukemic_seurat@meta.data %>% filter(project %in% c("T-LGLL", "Healthy")) %>% pull(barcode)
immune_seurat_tlgll_healthy <- subset(immune_seurat3, cells = cells.to.keep)
immune_seurat_tlgll_healthy <- immune_seurat_tlgll_healthy %>% preprocessSeurat(cells = colnames(immune_seurat_tlgll_healthy))

immune_tlgll_healthy_deg <- lapply(unique(immune_seurat_tlgll_healthy$cluster), getDEGbyClusterLGL, seurat_object = immune_seurat_tlgll_healthy) %>% rbindlist()

hallmark_up_lgll_v_healthy_up <- lapply(unique(immune_tlgll_healthy_deg$cluster), FUN = function(x) getHypergeometric(genes_df = subset(immune_tlgll_healthy_deg, avg_logFC > 0 & cluster == x), universe_df = universe_df, term_df = hallmark) %>% mutate(cluster = x)) %>% rbindlist()
hallmark_up_lgll_v_healthy_sigf_up <- hallmark_up_lgll_v_healthy_up %>% filter(p.adjust < 0.05) %>% mutate(ID = gsub("HALLMARK_", "", ID)) %>% mutate(ID = gsub("\\_", "\\ ", ID))

hallmark_dn_lgll_v_healthy_dn <- lapply(unique(immune_tlgll_healthy_deg$cluster), FUN = function(x) getHypergeometric(genes_df = subset(immune_tlgll_healthy_deg, avg_logFC < 0 & cluster == x), universe_df = universe_df, term_df = hallmark) %>% mutate(cluster = x)) %>% rbindlist()
hallmark_dn_lgll_v_healthy_sigf_dn <- hallmark_dn_lgll_v_healthy_dn %>% filter(p.adjust < 0.05) %>% mutate(ID = gsub("HALLMARK_", "", ID)) %>% mutate(ID = gsub("\\_", "\\ ", ID))


## b) cancer to healthy
cells.to.keep <- non_leukemic_seurat@meta.data %>% filter(!project %in% c("T-LGLL")) %>% pull(barcode)
immune_seurat_cancer_healthy <- subset(non_leukemic_seurat, cells = cells.to.keep)
immune_seurat_cancer_healthy <- immune_seurat_cancer_healthy %>% preprocessSeurat(cells = colnames(immune_seurat_cancer_healthy))

non_leukemic_lgll_v_cancer <- lapply(unique(immune_seurat_cancer_healthy$cluster), getDEGbyClusterLGL, seurat_object = immune_seurat_cancer_healthy) %>% rbindlist()

non_leukemic_lgll_v_cancer_up_hallmark <- lapply(unique(non_leukemic_lgll_v_cancer$cluster), FUN = function(x) getHypergeometric(genes_df = subset(non_leukemic_lgll_v_cancer, avg_logFC > 0 & cluster == x), universe_df = universe_df, term_df = hallmark) %>% mutate(cluster = x)) %>% rbindlist() %>% filter(p.adjust < 0.05)
non_leukemic_lgll_v_cancer_dn_hallmark <- lapply(unique(non_leukemic_lgll_v_cancer$cluster), FUN = function(x) getHypergeometric(genes_df = subset(non_leukemic_lgll_v_cancer, avg_logFC < 0 & cluster == x), universe_df = universe_df, term_df = hallmark) %>% mutate(cluster = x)) %>% rbindlist() %>% filter(p.adjust < 0.05)

non_leukemic_lgll_v_cancer_up_go <- lapply(unique(non_leukemic_lgll_v_cancer$cluster), FUN = function(x) getHypergeometric(genes_df = subset(non_leukemic_lgll_v_cancer, avg_logFC > 0 & cluster == x), universe_df = universe_df, term_df = go) %>% mutate(cluster = x)) %>% rbindlist() %>% filter(p.adjust < 0.05)
non_leukemic_lgll_v_cancer_dn_go <- lapply(unique(non_leukemic_lgll_v_cancer$cluster), FUN = function(x) getHypergeometric(genes_df = subset(non_leukemic_lgll_v_cancer, avg_logFC < 0 & cluster == x), universe_df = universe_df, term_df = go) %>% mutate(cluster = x)) %>% rbindlist() %>% filter(p.adjust < 0.05)
