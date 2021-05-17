

## Non-leukemic immune cells and T-LGLL vs 
## non-leukemic immunecells and hyperexpanded clonotypes

################################ Init cellphone input files

## Prepare for cellphonedb
cells.to.keep                  <- immune_seurat3@meta.data %>% filter(project %in% c("T-LGLL", "Healthy")) %>% pull(barcode)
lgl_healthy_from_immune_seurat <- subset(lgl_healthy_from_immune_seurat, cells = cells.to.keep)


## T-LGLL; take the immune-seurat and add the T-LGLL clones there
cells.to.keep          <- immune_seurat3@meta.data %>% filter(project %in% c("T-LGLL")) %>% pull(barcode)
lgl_from_immune_seurat <- subset(immune_seurat3, cells = cells.to.keep)
lgl_from_immune_seurat <- merge(lgl_from_immune_seurat, lgl_clones_seurat)
lgl_from_immune_seurat <- lgl_from_immune_seurat %>% preprocessSeurat(cells = colnames(lgl_from_immune_seurat))

lgl_from_immune_seurat$cellphone <- Idents(lgl_from_immune_seurat) %>% as.character()
lgl_from_immune_seurat$cellphone <- ifelse(!is.na(lgl_from_immune_seurat$clone_dipa), paste0(getId(lgl_from_immune_seurat$clone_dipa)), lgl_from_immune_seurat$cellphone)

table(lgl_from_immune_seurat$cellphone)

lgl_from_immune_seurat$cluster <- lgl_from_immune_seurat$cellphone
data.frame(lgl_from_immune_seurat@meta.data, cluster = lgl_from_immune_seurat$cellphone) %>% 
  initCellphonedb(seurat_object = lgl_from_immune_seurat, 
                  name = "tlgll_v_non_leukemic", 
                  sample_n = 50, 
                  folder = "results/manuscript/cellphone/")


## Healthy; take healthy from the immune-seruat and call out the hyperexpanded
cells.to.keep              <- lgl_healthy_from_immune_seurat@meta.data %>% filter(project %in% c("Healthy")) %>% pull(barcode)
healthy_from_immune_seurat <- subset(immune_seurat3, cells = cells.to.keep)
healthy_from_immune_seurat <- healthy_from_immune_seurat %>% preprocessSeurat(cells = colnames(healthy_from_immune_seurat))

# hyperexp_cells <- lgl_healthy_tcr_seurat3@meta.data %>% filter(project == "Healthy") %>% pull(barcode)
hyperexp_cells <- lgl_healthy_tcr_seurat3@meta.data %>% filter(project == "Healthy") %>% dplyr::select(barcode, new_clonotypes_id)
df <- healthy_from_immune_seurat@meta.data %>% left_join(hyperexp_cells, by = "barcode")
df$new_clonotypes_id.y <- gsub("_clonotype", " Cl ", df$new_clonotypes_id.y)
rownames(df) <- df$barcode
healthy_from_immune_seurat@meta.data <- df

healthy_from_immune_seurat$cellphone <- Idents(healthy_from_immune_seurat) %>% as.character()
healthy_from_immune_seurat$cellphone <- ifelse(!is.na(healthy_from_immune_seurat$new_clonotypes_id.y), healthy_from_immune_seurat$new_clonotypes_id.y, healthy_from_immune_seurat$cellphone)

healthy_from_immune_seurat$cluster <- healthy_from_immune_seurat$cellphone

dir.create("results/manuscript/cellphone/hyperexpanded_v_non_leukemic", showWarnings = F)
data.frame(healthy_from_immune_seurat@meta.data, cluster = healthy_from_immune_seurat$cellphone) %>% 
  initCellphonedb(seurat_object = healthy_from_immune_seurat, 
                  name = "hyperexpanded_v_non_leukemic", 
                  sample_n = 50, 
                  folder = "results/manuscript/cellphone/")



################################ Analyze

sigf_means_lgl     <- fread("results/manuscript/cellphone/out/tlgll_v_non_leukemic/significant_means.txt")
sigf_means_healthy <- fread("results/manuscript/cellphone/out/hyperexpanded_v_non_leukemic/significant_means.txt")

## How many interactions are there?
sigf_means_lgl_pairs <- getPairAmounts(sigf_means_lgl) %>% filterRedundantPairs()#  %>% mutate(cluster1 = getNewClusters(cluster1), cluster2 = getNewClusters(cluster2))
sigf_means_healthy_pairs <- getPairAmounts(sigf_means_healthy) %>% filterRedundantPairs()#  %>% mutate(cluster1 = getNewClusters(cluster1), cluster2 = getNewClusters(cluster2))
sigf_means <- rbind(sigf_means_lgl_pairs, sigf_means_healthy_pairs)

p <- cast(sigf_means, cluster1 ~ cluster2, sum)
rownames(p) <- colnames(p)[-1]
p <- p[,-1]
p[is.na(p)] <- 0

pdf("results/manuscript/cellphone/heatmap_healthy_lgl_clonotypes_tot2.pdf", width = 9, height = 5)
set.seed(123)
Heatmap(p,
        col=colorRamp2(c(0, 20, 40), c("dodgerblue4", "white", "darkred")),
        column_title_gp = gpar(fontsize = 8),
        row_title_gp = gpar(fontsize = 8),
        row_names_gp = gpar(fontsize = 9),
        column_names_gp = gpar(fontsize = 9),
        column_km = 4, border = T, 
        cluster_rows = T, 
        name = "#interactions")
dev.off()


