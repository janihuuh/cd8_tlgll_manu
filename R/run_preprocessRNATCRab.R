
message("Merge seurat-objects with TCRab-info...")

mergeTCRtoSeurat <- function(seurat_object, tcr_df){

  ## Merge with TCRab data with seurat-object metadata
  metadata_temp <- merge(seurat_object@meta.data, tcr_df, all.x = T, by.x = "barcode", by.y = "barcode_uniq")
  metadata_temp <- metadata_temp[match(colnames(seurat_object), metadata_temp$barcode), ]
  
  ## Add some meta data;
  
  ## 1) Major: over 10 cells in clonotype
  major_clonotypes               <- unique(subset(metadata_temp, metadata_temp$frequency > 10)$new_clonotypes_id)
  metadata_temp$major_clonotypes <- metadata_temp$new_clonotypes_id
  metadata_temp$major_clonotypes[!metadata_temp$new_clonotypes_id %in% major_clonotypes] <- "minor"
  
  
  ## 2) Expanded: over 2 cells in clonotype
  expanded_clonotypes <- unique(subset(metadata_temp, metadata_temp$frequency > 2)$new_clonotypes_id)
  metadata_temp$expanded_clonotypes <- metadata_temp$new_clonotypes_id
  metadata_temp$expanded_clonotypes[!metadata_temp$new_clonotypes_id %in% expanded_clonotypes] <- "unexpanded"
  
  
  ## Add metadata into Seurat object; make sure that the colnames match
  rownames(metadata_temp) <- metadata_temp$barcode
  colnames(seurat_object) == rownames(metadata_temp)
  
  
  seurat_object@meta.data <- metadata_temp
  return(seurat_object)
  
}

metadata <- lgl_seurat@meta.data[,1:52]
rownames(metadata) <- metadata$barcode
lgl_seurat@meta.data <- metadata

lgl_seurat <- removeTCRabData(lgl_seurat)
lgl_seurat <- mergeTCRtoSeurat(seurat_object = lgl_seurat, tcr_df = tot_barcode)

saveRDS(lgl_seurat, "results/lgl_seurat_new.rds")

tot_barcode <- fread("data/scRNAseq+TCRseq/preprocessed/tot_healthy_barcode.txt")
healthy_seurat <- removeTCRabData(healthy_seurat)
healthy_seurat <- mergeTCRtoSeurat(seurat_object = healthy_seurat, tcr_df = tot_barcode)


saveRDS(lgl_seurat, "results/lgl_seurat_new.rds")




tcr_df %>% filter(new_clonotypes_id %in% "FM805_clonotype1") %>% pull(barcode_uniq) %>% substr(start = 1, stop = 10) %>% table()
tcr_df %>% filter(new_clonotypes_id %in% "FM805_clonotype1") %>% View


df <- metadata_temp %>% filter(new_clonotypes_id %in% "FM805_clonotype1")
seurat_object$barcode.y <- NULL
barcodes1 <- substr(df$barcode.y, 1, nchar(df$barcode.y) - 2)
barcodes2 <- substr(df$barcode, nchar(df$barcode) - 15, nchar(df$barcode))

table(barcodes1 == barcodes2)
