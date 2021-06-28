

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



matchFile <- function(file, toMatch){

  message(extractFileName(file))
  # df <- fread(file) %>% mutate(name = extractFileName(file)) %>% mutate(key = paste(v,cdr3aa,j)) %>% dplyr::select(freq,count,name,key)
  # df <- merge(lgl_clones, df, by = "key")

  df <- fread(file) %>% mutate(name = extractFileName(file)) %>% dplyr::select(freq,count,name,cdr3aa) %>% group_by(cdr3aa,name) %>% summarise(count = sum(count), freq = sum(freq))
  df <- merge(toMatch, df, by = "cdr3aa")

}


getClusterPhenotypesLGLLvHealthyTCR2 <- function(clusters){

  clusters <- plyr::revalue(clusters, replace = c(

    "0"  = "0 CD8 T-LGLL",
    "1"  = "1 CD4 CM/naive" ,
    "2"  = "2 CD8 EM" ,
    "3"  = "3 CD4 EM" ,
    "4"  = "4 CD8 memory" ,
    "5"  = "5 CD8 T-LGLL" ,
    "6"  = "6 CD4" ,
    "7"  = "7 CD4 treg" ,
    "8"  = "8 CD4" ,
    "9"  = "9 CD8" ,
    "10" = "10 CD8 T-LGLL" ,
    "11" = "11 CD8 T-LGLL" ,
    "12" = "12 CD4 memory" ,
    "13" = "13 CD4" ,
    "14" = "14 CD8 T-LGLL" ))

  return(clusters)

}





vdjToGliph2 <- function(x){

  x %>% dplyr::select(cdr3aa, v, j, name, freq) %>%
    dplyr::rename("CDR3b"   = "cdr3aa",
                  "TRBV"    = "v",
                  "TRBJ"    = "j",
                  "Patient" =  name,
                  "Counts"  = freq) %>% mutate(CDR3a = "NA") %>% dplyr::select(CDR3b, TRBV, TRBJ, CDR3a, Patient, Counts) %>% dplyr::rename("subject:condition" = "Patient",	"count" = "Counts")

}


writeGliph2Commands <- function(filename, output_folder){

  df <- paste0(paste0("out_prefix=", gsub(".txt", "", extractFileName(filename))), "\n",
               paste0("cdr3_file=", filename), "\n",
               "hla_file=/Users/janihuuh/Dropbox/applications/gliph2.2/data/aa_test_hla.txt", "\n",

               "refer_file=/Users/janihuuh/Dropbox/applications/gliph2.2/human_v2.0/ref_CD8_v2.0.txt", "\n",
               "v_usage_freq_file=/Users/janihuuh/Dropbox/applications/gliph2.2/human_v2.0/ref_V_CD8_v2.0.txt", "\n",
               "cdr3_length_freq_file=/Users/janihuuh/Dropbox/applications/gliph2.2/human_v2.0/ref_L_CD8_v2.0.txt", "\n",
               "local_min_pvalue=0.001", "\n",
               "p_depth=1000", "\n",
               "global_convergence_cutoff=1", "\n",
               "simulation_depth=1000", "\n",
               "kmer_min_depth=5", "\n",
               "local_min_OVE=10", "\n",
               "algorithm=GLIPH2", "\n",
               "all_aa_interchangeable=1")

  fwrite(as.data.frame(df), paste0(output_folder, gsub(".txt", ".sh", extractFileName(filename))), sep = "\t", quote = F, row.names = F, col.names = F)}







getClusterPhenotypesLGLLvHealthy_ver1 <- function(clusters){

  clusters <- plyr::revalue(clusters, replace = c(

    "0"  = "0 CD8 effector",
    "1"  = "1 CD4" ,
    "2"  = "2 CD8 T-LGLL" ,
    "3"  = "3 Monocytes CD16-" ,
    "4"  = "4 T-cells" ,
    "5"  = "5 NK-cells" ,
    "6"  = "6 B-cells" ,
    "7"  = "7 Monocytes CD14+" ,
    "8"  = "8 CD8 T-LGLL" ,
    "9"  = "9 Other" ,
    "10" = "10 CD8" ,
    "11" = "11 Monocytes CD16+" ,
    "12" = "12 pDC" ))

  return(clusters)

}



getClusterPhenotypesCompare <- function(clusters){

  clusters <- plyr::revalue(clusters, replace = c(

    "0"  = "0 CD8 EM",
    "1"  = "1 Monocyte CD16-" ,
    "2"  = "2 CD4 EM" ,
    "3"  = "3 NK CD56dim" ,
    "4"  = "4 CD4 naive" ,
    "5"  = "5 T batch" ,
    "6"  = "6 T batch" ,
    "7"  = "7 Monocyte CD16+" ,
    "8"  = "8 CD4 CM" ,
    "9"  = "9 B memory" ,
    "10" = "10 Monocytes" ,
    "11" = "11 B naive" ,
    "12" = "12 CD8 naive" ,
    "13" = "13 T-cell naive/CM" ,
    "14" = "14 CD8 effector" ,
    "15" = "15 CD4 EM" ,
    "16" = "16 cDC" ,
    "17" = "17 CD8 effector" ,
    "18" = "18 B naive" ,
    "19" = "19 Neutrophils" ,
    "20" = "20 pDC" ,
    "21" = "21 NK cycling" ,
    "22" = "22 B plasma"))

  return(clusters)

}





gliphExtractClusteredTCRs <- function(df){

  i <- 1
  df_list <- NULL

  for(i in 1:nrow(df)){

    x = df[i,]
    df_list[[i]] <- strsplit(noquote(x$V3), split = "\\ ")[[1]]

  }

  df <- df_list %>% do.call(what = "c")
  return(df)

}


gliphDfToEdgelist <- function(df){

  i <- 1
  df_list <- NULL

  for(i in 1:nrow(df)){

    x = df[i,]
    df_list[[i]] <- data.frame(V1 = gsub(x = x$V2, pattern = "CRG\\-", replacement = ""),
                               V2 = strsplit(noquote(x$V3), split = "\\ ")[[1]])

  }

  df <- df_list %>% rbindlist()
  return(df)

}



plotGliphRemix <- function(el, color){

  el <- el %>% as.matrix()
  el <- graph_from_edgelist(el)
  g_am <- simplify(el, remove.multiple = F, remove.loops = T)


  ## Plot
  V(g_am)$color <- ifelse(V(g_am) %>% names() %in% lgl_annotation$cdr3aa == T, "darkblue", color)
  V(g_am)$frame.color <- ifelse(V(g_am) %>% names() %in% lgl_annotation$cdr3aa == T, "darkblue", color)
  V(g_am)$size <- ifelse(V(g_am) %>% names() %in% lgl_annotation$cdr3aa == T, 5, 1.5)
  V(g_am)$label <- ""

  E(g_am)$arrow.mode <- 0
  E(g_am)$color <- "gray80"
  E(g_am)$size <- 1

  layout = layout_with_graphopt(g_am)

  par(mar=c(0,0,0,0)+.1)
  plot.igraph(g_am, layout = layout, axes = F, xlim = c(-0.5,0.5), add = F) #, ylim = c(-0.5,0.5))

}






getDoublets <- function(seurat_object){

  require(scds)

  # Annotate doublet using co-expression based doublet scoring:
  sce_object <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = seurat_object@assays$RNA@counts), colData = seurat_object@meta.data)
  sce_object <- cxds(sce_object)
  sce_object <- bcds(sce_object)
  sce_object <- cxds_bcds_hybrid(sce_object)

  ## Add into Seurat
  seurat_object$cxds_doublet_score   <- SingleCellExperiment::colData(sce_object)$cxds_score
  seurat_object$bcds_doublet_score   <- SingleCellExperiment::colData(sce_object)$bcds_score

  seurat_object$hybrid_doublet_score <- SingleCellExperiment::colData(sce_object)$hybrid_score
  seurat_object$cxds_doublet_score_norm <- c(SingleCellExperiment::colData(sce_object)$cxds_score - min(SingleCellExperiment::colData(sce_object)$cxds_score)) / max(SingleCellExperiment::colData(sce_object)$cxds_score)
  seurat_object$bcds_doublet_score_norm <- c(SingleCellExperiment::colData(sce_object)$bcds_score - min(SingleCellExperiment::colData(sce_object)$bcds_score)) / max(SingleCellExperiment::colData(sce_object)$bcds_score)
  return(seurat_object)

}



getClustering <- function(seurat_object){

  ## Clustering
  res           <- c(seq(0.1, 1, 0.1), seq(1.2, 2, 0.2), 2.5, 3)
  seurat_object <- FindNeighbors(seurat_object, reduction = "pca", dims = c(1:ncol(seurat_object@reductions$pca@cell.embeddings)))
  seurat_object <- FindClusters(object = seurat_object, resolution = res, verbose = F)
  return(seurat_object)

}


plotQC <- function(seurat_object, folder){

  # min_mito     <- 0
  # max_mito     <- 15
  #
  # min_ribo     <- 5
  # max_ribo     <- 50
  #
  # min_features <- 300
  # max_features <- 6e3
  #
  # min_counts   <- 1e3
  # max_counts   <- 45e3

  min_mito     <- 0
  max_mito     <- 15

  min_ribo     <- 5
  max_ribo     <- 50

  min_features <- 300
  max_features <- 5e3

  min_counts   <- 1e3
  max_counts   <- 30e3

  qc_df <- seurat_object@meta.data %>% as.data.frame()
  qc_df$cluster = Idents(seurat_object)

  plotQcViolin(qc_df, var_to_plot = "nFeature_RNA", grouping = "cluster", min = min_features, max = max_features) + theme_classic(base_size = 17) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(paste0(folder, "violin_nFeature_RNA.png"), width = 6, height = 4)

  plotQcViolin(qc_df, var_to_plot = "nCount_RNA", grouping = "cluster", min = min_counts, max = max_counts) + theme_classic(base_size = 17) + scale_y_log10() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(paste0(folder, "violin_nCount_RNA.png"), width = 6, height = 4)

  plotQcViolin(qc_df, var_to_plot = "percent.mt", grouping = "cluster", min = min_mito, max = max_mito) + theme_classic(base_size = 17) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(paste0(folder, "violin_percent_mt.png"), width = 6, height = 4)

  plotQcViolin(qc_df, var_to_plot = "percent.ribo", grouping = "cluster", min = min_ribo, max = max_ribo) + theme_classic(base_size = 17) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(paste0(folder, "/violin_percent_ribo.png"), width = 6, height = 4)

  ## Scatter plots
  p <- qc_df %>%
    ggplot(aes(nCount_RNA, nFeature_RNA, color = orig.ident)) + geom_point(size = 0.3, alpha = 0.5) + scale_x_log10() + scale_y_log10() +
    geom_vline(xintercept = min_counts, linetype = "dotted") +
    geom_vline(xintercept = max_counts, linetype = "dotted") +
    geom_hline(yintercept = min_features, linetype = "dotted") +
    geom_hline(yintercept = max_features, linetype = "dotted") + theme(legend.position = "none")
  ggsave(plot = p, paste0(folder, "/scatter_counts_vs_genes.png"), width = 5, height = 4)

  p <- qc_df %>%
    ggplot(aes(percent.mt, nCount_RNA, color = orig.ident)) + geom_point(size = 0.3, alpha = 0.5) + add_guide +
    scale_y_log10() +
    geom_hline(yintercept = min_counts, linetype = "dotted") +
    geom_hline(yintercept = max_counts, linetype = "dotted") +
    geom_vline(xintercept = max_mito, linetype = "dotted") + theme(legend.position = "none")
  ggsave(plot = p, paste0(folder,"/scatter_mito_vs_counts.png"), width = 5, height = 4)

  p <- qc_df %>%
    ggplot(aes(percent.mt, nFeature_RNA, color = orig.ident)) + geom_point(size = 0.3, alpha = 0.5) + add_guide +
    scale_y_log10() +
    geom_hline(yintercept = min_features, linetype = "dotted") +
    geom_hline(yintercept = max_features, linetype = "dotted") +
    geom_vline(xintercept = max_mito, linetype = "dotted") + theme(legend.position = "none")
  ggsave(plot = p, paste0(folder,"/scatter_mito_vs_genes.png"), width = 5, height = 4)

  p <- qc_df %>%
    ggplot(aes(percent.ribo, percent.mt, color = orig.ident)) + geom_point(size = 0.3, alpha = 0.5) + add_guide +
    geom_hline(yintercept = min_mito, linetype = "dotted") +
    geom_hline(yintercept = max_mito, linetype = "dotted") +
    geom_vline(xintercept = min_ribo, linetype = "dotted") +
    geom_vline(xintercept = max_ribo, linetype = "dotted") + theme(legend.position = "none")
  ggsave(plot = p, paste0(folder, "/scatter_ribo_vs_mito.png"), width = 5, height = 4)

}



plotQC_manu <- function(seurat_object, folder){

  plotQcViolin2 <- function(viz_df, var_to_plot, grouping, min, max){

    viz_df_temp <- viz_df %>% dplyr::select(var_to_plot)

    label_df_min <- ifelse(viz_df_temp > min, "above", "below") %>% table
    label_df_max <- ifelse(viz_df_temp < max, "above", "below") %>% table

    ggplot(data = viz_df, aes_string(x = grouping, y = var_to_plot, fill = grouping)) +
      geom_violin(alpha = 0.5) +
      # geom_boxplot(width = 0.1, alpha = 0.8, outlier.shape = NA) +

      # geom_hline(yintercept = min, linetype = "dotted") +
      # geom_hline(yintercept = max, linetype = "dotted") +

      # annotate(geom = "text", x = 2.5, y = min, label = paste("Below the line:\n", label_df_min[2]), fontface = "italic") +
      # annotate(geom = "text", x = 2.5, y = max, label = paste("Above the line:\n", label_df_max[2]), fontface = "italic") +

      labs(x = "", title = var_to_plot) + theme(legend.position = "none")

  }

  min_mito     <- 0
  max_mito     <- 15

  min_ribo     <- 0
  max_ribo     <- 50

  min_features <- 300
  max_features <- 5e3

  min_counts   <- 1e3
  max_counts   <- 30e3

  qc_df <- seurat_object@meta.data %>% as.data.frame()
  qc_df$cluster = Idents(seurat_object)

  plotQcViolin2(qc_df, var_to_plot = "nFeature_RNA", grouping = "cluster", min = 0, max = max_features) + theme_classic(base_size = 17) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(legend.position = "none")
  ggsave(paste0(folder, "violin_nFeature_RNA.png"), width = 6, height = 4)

  plotQcViolin2(qc_df, var_to_plot = "nCount_RNA", grouping = "cluster", min = 0, max = max_counts) + theme_classic(base_size = 17) + scale_y_log10() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(legend.position = "none")
  ggsave(paste0(folder, "violin_nCount_RNA.png"), width = 6, height = 4)

  plotQcViolin2(qc_df, var_to_plot = "percent.mt", grouping = "cluster", min = 0, max = max_mito) + theme_classic(base_size = 17) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(legend.position = "none")
  ggsave(paste0(folder, "violin_percent_mt.png"), width = 6, height = 4)

  plotQcViolin2(qc_df, var_to_plot = "percent.ribo", grouping = "cluster", min = 0, max = max_ribo) + theme_classic(base_size = 17) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(legend.position = "none")
  ggsave(paste0(folder, "/violin_percent_ribo.png"), width = 6, height = 4)

}







getId <- function(vec){

  vec <- gsub("FM805", "LGLL_1", vec)
  vec <- gsub("FM1656", "LGLL_2", vec)
  vec <- gsub("FM1382", "LGLL_3", vec)
  vec <- gsub("FM1338", "LGLL_4", vec)
  vec <- gsub("FM1890", "LGLL_5", vec)
  vec <- gsub("FM1189", "LGLL_6", vec)
  vec <- gsub("FM2481", "LGLL_7", vec)
  vec <- gsub("FM2477", "LGLL_8", vec)
  vec <- gsub("FM2311", "LGLL_9", vec)
  return(vec)
}

getQCSzaboT <- function(seurat_object){

  ###################

  min_mito     <- 0
  max_mito     <- 15

  min_ribo     <- 10
  max_ribo     <- 50

  min_features <- 250
  max_features <- 3e3

  min_counts   <- 1000
  max_counts   <- 10e3

  cycle.genes  <- c("ANLN", "ASPM","BIRC5","CCNA2","CCNB1","CCNB2","CCND1","CD63","CDC20","CDCA8","CDKN3","CENPE","CENPF",
                    "CEP55","CKAP2L","DLGAP5","FOXM1","GTSE1","H2AFZ","HIST1H1B", "HIST1H1C", "HIST1H1D", "HIST1H1E", "HIST1H2AJ",
                    "HIST1H4C", "HJURP", "HMGB1", "HMGB2", "HMMR", "KIF11", "KIF14", "KIF15", "KIF2C", "LMNA",
                    "MCM3", "MKI67", "NCAPG", "NUSAP1", "PCNA", "PLK1", "PRC1", "RRM2", "SMC4", "STMN1", "TK1", "TOP2A", "TPX2", "TUBA1B",
                    "TUBB", "TYMS", "UBE2C")

  ###################

  seurat_object  <- PercentageFeatureSet(seurat_object, pattern = "^MT-", col.name = "percent.mt")
  seurat_object  <- PercentageFeatureSet(seurat_object, pattern = "^RP", col.name = "percent.ribo")
  # seurat_object  <- PercentageFeatureSet(seurat_object, features = cycle.genes, col.name = "percent.cycle")
  seurat_object@meta.data$barcode <- colnames(seurat_object)

  ## In total, we remove with the following conditions:
  qc_df <- seurat_object@meta.data %>% as.data.frame()

  # VlnPlot(seurat_object, group.by = "Patient_ID", features = "percent.mt") # + geom_hline(yintercept = 4500) + geom_hline(yintercept = 100)
  # VlnPlot(seurat_object, group.by = "Patient_ID", features = "percent.ribo") # + geom_hline(yintercept = 4500) + geom_hline(yintercept = 100)
  # VlnPlot(seurat_object, group.by = "Patient_ID", features = "nFeature_RNA") # + geom_hline(yintercept = 4500) + geom_hline(yintercept = 100)
  # VlnPlot(seurat_object, group.by = "Patient_ID", features = "nCount_RNA") + geom_hline(yintercept = 30e3) + geom_hline(yintercept = 5e2)


  percent_mito_outlier <- qc_df %>% dplyr::filter(percent.mt   > max_mito     | percent.mt   < min_mito)     %>% pull(barcode) %>% as.character()
  percent_ribo_outlier <- qc_df %>% dplyr::filter(percent.ribo > max_ribo     | percent.ribo < min_ribo)     %>% pull(barcode) %>% as.character()
  features_outlier     <- qc_df %>% dplyr::filter(nFeature_RNA < min_features | nFeature_RNA > max_features) %>% pull(barcode) %>% as.character()
  umis_outlier         <- qc_df %>% dplyr::filter(nCount_RNA   > max_counts   | nCount_RNA   < min_counts)   %>% pull(barcode) %>% as.character()

  outlier_cells        <- c(percent_mito_outlier,
                            percent_ribo_outlier,
                            features_outlier,
                            umis_outlier)

  reason               <- c(rep("percent_mito_outlier", length(percent_mito_outlier)),
                            rep("percent_ribo_outlier", length(percent_ribo_outlier)),
                            rep("features_outlier",     length(features_outlier)),
                            rep("umis_outlier",         length(umis_outlier)))

  outlier_df <- data.frame(barcode = outlier_cells, reason = reason) %>% dplyr::mutate(from = extractName(barcode)) #, 1, 10))

  ## Remove the cells from Seurat-object and save a new seurat-object
  cells.to.use  <- colnames(seurat_object)[!colnames(seurat_object) %in% outlier_df$barcode]
  seurat_object <- subset(seurat_object, cells = cells.to.use)
  return(seurat_object)

}

getScviInput <- function(seurat_object, folder){

  dir.create(folder, showWarnings = F)
  genes_to_keep <- list()
  i <- 1
  Idents(seurat_object) <- seurat_object$orig.ident

  for(patient in unique(seurat_object$orig.ident)){

    message(patient)
    seurat_temp <- subset(seurat_object, idents = patient)
    counts_temp <- seurat_temp@assays$RNA@counts %>% as.data.frame

    genes_to_keep[[i]] <- rownames(counts_temp)
    i <- i + 1

    counts_temp <- seurat_object@assays$RNA@counts[ ,seurat_object$orig.ident == patient] %>% as.data.frame
    counts_temp <- counts_temp[!rownames(counts_temp) %in% clonality_genes, ]
    data.table::fwrite(counts_temp, paste0(folder, patient, ".csv"), sep = ",", quote = F, row.names = T, col.names = T)

  }
}



getQC <- function(seurat_object){

  ###################

  min_mito     <- 0
  max_mito     <- 15

  min_ribo     <- 10
  max_ribo     <- 50

  min_features <- 250
  max_features <- 4500

  min_counts   <- 1000
  max_counts   <- 20e3

  cycle.genes  <- c("ANLN", "ASPM","BIRC5","CCNA2","CCNB1","CCNB2","CCND1","CD63","CDC20","CDCA8","CDKN3","CENPE","CENPF",
                    "CEP55","CKAP2L","DLGAP5","FOXM1","GTSE1","H2AFZ","HIST1H1B", "HIST1H1C", "HIST1H1D", "HIST1H1E", "HIST1H2AJ",
                    "HIST1H4C", "HJURP", "HMGB1", "HMGB2", "HMMR", "KIF11", "KIF14", "KIF15", "KIF2C", "LMNA",
                    "MCM3", "MKI67", "NCAPG", "NUSAP1", "PCNA", "PLK1", "PRC1", "RRM2", "SMC4", "STMN1", "TK1", "TOP2A", "TPX2", "TUBA1B",
                    "TUBB", "TYMS", "UBE2C")

  ###################

  seurat_object  <- PercentageFeatureSet(seurat_object, pattern = "^MT-", col.name = "percent.mt")
  seurat_object  <- PercentageFeatureSet(seurat_object, pattern = "^RP", col.name = "percent.ribo")
  # seurat_object  <- PercentageFeatureSet(seurat_object, features = cycle.genes, col.name = "percent.cycle")
  seurat_object@meta.data$barcode <- colnames(seurat_object)

  ## In total, we remove with the following conditions:
  qc_df <- seurat_object@meta.data %>% as.data.frame()

  # VlnPlot(seurat_object, group.by = "Patient_ID", features = "percent.mt") # + geom_hline(yintercept = 4500) + geom_hline(yintercept = 100)
  # VlnPlot(seurat_object, group.by = "Patient_ID", features = "percent.ribo") # + geom_hline(yintercept = 4500) + geom_hline(yintercept = 100)
  # VlnPlot(seurat_object, group.by = "Patient_ID", features = "nFeature_RNA") # + geom_hline(yintercept = 4500) + geom_hline(yintercept = 100)
  # VlnPlot(seurat_object, group.by = "Patient_ID", features = "nCount_RNA") + geom_hline(yintercept = 30e3) + geom_hline(yintercept = 5e2)


  percent_mito_outlier <- qc_df %>% dplyr::filter(percent.mt   > max_mito     | percent.mt   < min_mito)     %>% pull(barcode) %>% as.character()
  percent_ribo_outlier <- qc_df %>% dplyr::filter(percent.ribo > max_ribo     | percent.ribo < min_ribo)     %>% pull(barcode) %>% as.character()
  features_outlier     <- qc_df %>% dplyr::filter(nFeature_RNA < min_features | nFeature_RNA > max_features) %>% pull(barcode) %>% as.character()
  umis_outlier         <- qc_df %>% dplyr::filter(nCount_RNA   > max_counts   | nCount_RNA   < min_counts)   %>% pull(barcode) %>% as.character()

  outlier_cells        <- c(percent_mito_outlier,
                            percent_ribo_outlier,
                            features_outlier,
                            umis_outlier)

  reason               <- c(rep("percent_mito_outlier", length(percent_mito_outlier)),
                            rep("percent_ribo_outlier", length(percent_ribo_outlier)),
                            rep("features_outlier",     length(features_outlier)),
                            rep("umis_outlier",         length(umis_outlier)))

  outlier_df <- data.frame(barcode = outlier_cells, reason = reason) %>% dplyr::mutate(from = extractName(barcode)) #, 1, 10))

  ## Remove the cells from Seurat-object and save a new seurat-object
  cells.to.use  <- colnames(seurat_object)[!colnames(seurat_object) %in% outlier_df$barcode]
  seurat_object <- subset(seurat_object, cells = cells.to.use)
  return(seurat_object)

}

getScviInput <- function(seurat_object, folder){

  dir.create(folder, showWarnings = F)
  genes_to_keep <- list()
  i <- 1
  Idents(seurat_object) <- seurat_object$orig.ident

  for(patient in unique(seurat_object$orig.ident)){

    message(patient)
    seurat_temp <- subset(seurat_object, idents = patient)
    counts_temp <- seurat_temp@assays$RNA@counts %>% as.data.frame

    genes_to_keep[[i]] <- rownames(counts_temp)
    i <- i + 1

    counts_temp <- seurat_object@assays$RNA@counts[ ,seurat_object$orig.ident == patient] %>% as.data.frame
    counts_temp <- counts_temp[!rownames(counts_temp) %in% clonality_genes, ]
    data.table::fwrite(counts_temp, paste0(folder, patient, ".csv"), sep = ",", quote = F, row.names = T, col.names = T)

  }
}



putLatentsSeurat <- function(seurat_object, latent){

  latent_umap <- uwot::umap(latent) %>% as.data.frame() %>% dplyr::rename(UMAP1=V1, UMAP2=V2)

  latent      <- as.matrix(latent)
  latent_umap <- as.matrix(latent_umap)

  rownames(latent)      <- colnames(seurat_object)
  rownames(latent_umap) <- colnames(seurat_object)

  latent_dim_red            <- CreateDimReducObject(key = "latent", embeddings = as.matrix(x = latent))
  latent_umap_dim_red       <- CreateDimReducObject(key = "latent_umap", embeddings = as.matrix(x = latent_umap))

  seurat_object[['latent']]      <- latent_dim_red
  seurat_object[['latent_umap']] <- latent_umap_dim_red
  return(seurat_object)
}


vdjToGliph <- function(vdj_df){

  ## Write gliph files to the vdj files
  # @ param
  # input: df from vdj

  df <- vdj_df %>% dplyr::select(cdr3aa, v, j, name, count) %>%
    dplyr::rename("CDR3b"   = "cdr3aa",
                  "TRBV"    = "v",
                  "TRBJ"    = "j",
                  "Patient" =  name,
                  "Counts"  = count)

  return(df)

}


hpca.se   <- SingleR::HumanPrimaryCellAtlasData()
blueprint <- SingleR::BlueprintEncodeData()

getSingler <- function(seurat_object, cluster = NULL, method = NULL, sample = NULL){

  if(!is.null(sample)){

    set.seed(123)
    seurat_object <- subset(seurat_object, cells = colnames(seurat_object)[sample(1:ncol(seurat_object), sample)])

  }

  sce       <- as.SingleCellExperiment(seurat_object)

  ## Predictions
  if(is.null(method)){
    pred.hca <- SingleR::SingleR(test = sce, ref = hpca.se, assay.type.test = 1,   labels = hpca.se$label.fine)
    pred.blu <- SingleR::SingleR(test = sce, ref = blueprint, assay.type.test = 1, labels = blueprint$label.fine)

    if(is.null(sample)){
      seurat_object$singler_hpca_pred      <- pred.hca$first.labels
      seurat_object$singler_blueprint_pred <- pred.hca$first.labels
      return(seurat_object)
    }

    else{
      df <- data.frame(barcode = rownames(pred.hca), cluster = seurat_object$cluster, singler_hpca_pred = pred.hca$labels, singler_blueprint_pred = pred.blu$labels)
      return(df)
    }

  }


  if(method == "cluster"){
    if(is.null(cluster)){
      cluster=seurat_object$cluster
    }
    pred.hca <- SingleR::SingleR(test = sce, ref = hpca.se, assay.type.test = 1,   labels = hpca.se$label.fine, method = "cluster", clusters = cluster)
    pred.blu <- SingleR::SingleR(test = sce, ref = blueprint, assay.type.test = 1, labels = blueprint$label.fine, method = "cluster", clusters = cluster)
    df <- data.frame(cluster = rownames(pred.hca), singler_hpca_pred = pred.hca$labels, singler_blueprint_pred = pred.blu$labels)
    return(df)
  }
}


plotClustering <- function(seurat_object){

  res <- c(seq(0.1, 2, 0.1), 2.5, 3)
  clustering_columns <- grep("res", colnames(seurat_object@meta.data), value = T)
  clustering_columns <- clustering_columns[order(substr(clustering_columns, 13, nchar(clustering_columns)) %>% as.numeric())]

  q <- NULL; i <- 1

  for(clustering_column in clustering_columns){
    q[[i]] <- seurat_object@meta.data[,clustering_column] %>% levels %>% length
    i <- i + 1
  }

  data.frame(resolution = res, nClusters = do.call(q, what="c")) %>%
    ggplot(aes((resolution),nClusters), label = nClusters) + geom_point(shape = 21) + theme_bw()

}


## Run seurat
preprocessSeurat <- function(orig_object, cells.to.use){

  ## Subset object
  object <- subset(orig_object, cells = cells.to.use)

  # orig_object@meta.data$barcode
  temp_meta <- orig_object@meta.data[as.character(orig_object@meta.data$barcode) %in% cells.to.use, ]
  temp_meta <- temp_meta[match(colnames(object), temp_meta$barcode), ]
  temp_meta$barcode == colnames(object)
  object@meta.data <- temp_meta

  ## Normalize and find HVGs
  object  <- NormalizeData(object, normalization.method = "LogNormalize", scale.factor = 10000)
  object  <- FindVariableFeatures(object, selection.method = "vst", nfeatures = 2000, clip.max = 10)

  ## Remove clonality genes
  hvg     <- VariableFeatures(object)
  too_hvg <- HVFInfo(object = object) %>% add_rownames(var = "gene") %>% filter(variance.standardized > 10) %>% pull("gene") %>% as.character()
  hvg     <- hvg[!hvg %in% too_hvg]
  hvg     <- hvg[!hvg %in% clonality_genes]
  hvg     <- hvg[!hvg %in% unwanted_genes]

  VariableFeatures(object) <- hvg

  ## Scale data
  object <- ScaleData(object, features = hvg)

  ## PCA data
  object <- RunPCA(object, features = hvg, npcs = 50)
  nPCs   <- sum(object[["pca"]]@stdev > 2)
  print(paste("nPCs:", nPCs))

  ## RunUMAP does not work
  object <- RunUMAP(object, dims = 1:nPCs, learning.rate = 1)
  return(object)

}

getDEGbyClusterLGL <- function(seurat_object, cluster, min_cells = 50){

  message(paste0("===== ", cluster, " ====="))

  ## If under 50 cells to begin with
  if(table(Idents(seurat_object)) %>% as.data.frame() %>% filter(Var1 == cluster) %>% pull(Freq) <= min_cells) return(NULL)

  ## Subet to only cluster
  seurat_cluster         <- subset(seurat_object, ident = cluster)
  Idents(seurat_cluster) <- seurat_cluster$lgl

  ## Calculate DEG only if there's at least 5 cells per time point
  n_df <- table(Idents(seurat_cluster)) %>% as.data.frame()

  n1 <- n_df %>% filter(Var1 == "T-LGLL") %>% pull(Freq) >= min_cells
  n2 <- n_df %>% filter(Var1 == "else")       %>% pull(Freq) >= min_cells

  if(length(n1) == 0) n1 <- FALSE
  if(length(n2) == 0) n2 <- FALSE

  if(n1 & n2){
    cluster_markers_2v1 <- FindMarkers(object = seurat_cluster, ident.1 = "else",  ident.2 = "T-LGLL",  only.pos = F, min.pct = 0.05, logfc.threshold = 0.1, return.thresh = 1, do.print = T, test.use = "t", max.cells.per.ident = 1e3) %>% add_rownames(var = "gene") %>% mutate(timepoint = "LGLLvHealthy") %>% mutate(p_val_adj = p.adjust(p_val, method = "BH"))
    df <- cluster_markers_2v1
    if(!is.null(df)) df <- df %>% filter(p_val_adj < 0.05) %>% mutate(cluster = cluster, direction = ifelse(avg_logFC > 0, "up", "down"))
    return(df)
  }
}

getNewClusters <- function(clusters){clusters %>% extractClusterNumber() %>% as.numeric() %>% as.factor() }

removeTCRabData <- function(seurat_object){

  seurat_object@meta.data$cdr3s_nt <- NULL
  seurat_object@meta.data$tra_cdr3s_nt <- NULL
  seurat_object@meta.data$trb_cdr3s_nt <- NULL
  seurat_object@meta.data$cdr3s_aa <- NULL
  seurat_object@meta.data$tra_cdr3s_aa <- NULL
  seurat_object@meta.data$trb_cdr3s_aa <- NULL
  seurat_object@meta.data$clonotype_id <- NULL
  seurat_object@meta.data$cdr3s_nt               <- NULL
  seurat_object@meta.data$tra_cdr3s_nt <- NULL
  seurat_object@meta.data$trb_cdr3s_nt <- NULL
  seurat_object@meta.data$cdr3s_aa <- NULL
  seurat_object@meta.data$tra_cdr3s_aa <- NULL
  seurat_object@meta.data$trb_cdr3s_aa <- NULL
  seurat_object@meta.data$clonotype_id <- NULL
  seurat_object@meta.data$frequency <- NULL
  seurat_object@meta.data$proportion <- NULL
  seurat_object@meta.data$cdr3s_aa.1 <- NULL
  seurat_object@meta.data$cdr3s_nt.1 <- NULL
  seurat_object@meta.data$chain_tra <- NULL
  seurat_object@meta.data$v_tra <- NULL
  seurat_object@meta.data$d_tra <- NULL
  seurat_object@meta.data$j_tra <- NULL
  seurat_object@meta.data$cdr3s_nt_freq_tra <- NULL
  seurat_object@meta.data$cdr3s_aa_freq_tra <- NULL
  seurat_object@meta.data$v_freq_tra <- NULL
  seurat_object@meta.data$d_freq_tra <- NULL
  seurat_object@meta.data$j_freq_tra <- NULL
  seurat_object@meta.data$chain_trb             <- NULL
  seurat_object@meta.data$v_trb <- NULL
  seurat_object@meta.data$d_trb <- NULL
  seurat_object@meta.data$j_trb <- NULL
  seurat_object@meta.data$cdr3s_nt_freq_trb <- NULL
  seurat_object@meta.data$cdr3s_aa_freq_trb <- NULL
  seurat_object@meta.data$v_freq_trb <- NULL
  seurat_object@meta.data$d_freq_trb <- NULL
  seurat_object@meta.data$j_freq_trb <- NULL
  seurat_object@meta.data$tcr_type <- NULL
  seurat_object@meta.data$patient <- NULL
  seurat_object@meta.data$new_clonotypes_id  <- NULL

  return(seurat_object)

}

getLatentClustering <- function(seurat_object){

  ## Clustering
  res <- c(seq(0.1, 2, 0.1), 2.5, 3)
  seurat_object <- FindNeighbors(seurat_object, reduction = "latent", dims = c(1:ncol(seurat_object@reductions$latent@cell.embeddings)))
  seurat_object <- FindClusters(object = seurat_object, resolution = res, verbose = F)
  return(seurat_object)

}

plotSlingshot <- function(slingshot_object, reducedDim){

  colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)

  cololors                <- slingshot_object$orig.clusters %>% extractClusterNumber() %>% as.numeric()
  cololors[cololors == 0] <- getPalette5(8)[1]
  cololors[cololors == 1] <- getPalette5(8)[2]
  cololors[cololors == 2] <- getPalette5(8)[3]
  cololors[cololors == 3] <- getPalette5(8)[4]
  cololors[cololors == 4] <- getPalette5(8)[5]
  cololors[cololors == 5] <- getPalette5(8)[6]
  cololors[cololors == 6] <- getPalette5(8)[7]
  cololors[cololors == 7] <- getPalette5(8)[8]


  curve_cols <- c("black", "darkred", "darkblue", "darkolivegreen4")

  colors2 <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(9)

  sling_curve <- SlingshotDataSet(slingshot_object)

  par(mfrow = c(1,2))
  plot(reducedDims(slingshot_object)[[reducedDim]], col = cololors, pch = 16, asp = 1)
  for(i in 1:2){
    lines(sling_curve@curves[[i]], lwd = 4, col = curve_cols[i])
  }

  plot(reducedDims(slingshot_object)[[reducedDim]], col = colors[cut(slingshot_object$slingPseudotime_1,breaks=100)], pch=16, asp = 1)
  # lines(SlingshotDataSet(slingshot_object), lwd = 2, col = "black")
  for(i in 1:2){
    lines(sling_curve@curves[[i]], lwd = 4, col = curve_cols[i])
  }

  par(mfrow = c(1,1))


}

runSlingshot <- function(sce_object, cluster_with_earliset_timepoint, reducedDim){

  require(slingshot); require(SummarizedExperiment)
  sce <- slingshot(data = sce_object, clusterLabels = 'orig.clusters', reducedDim = reducedDim, start.clus = cluster_with_earliset_timepoint)

}

getNewClusters <- function(clusters){
  clusters %>% extractClusterNumber() %>% as.numeric() %>% as.factor() %>% getClusterPhenotypesCml()
}

extractCoarsePhenotype <- function(strs){

  p <- NULL
  i <- 1
  for(str1 in strs){
    p[[i]] <- strsplit(str1, "[ ]")[[1]][2]
    i <- i + 1
  }

  return(p)

}

facets_nice <- theme(strip.background = element_rect(fill="grey96"), strip.text = element_text(colour = 'black'))

reorderClusters <- function(cluster_vec){

  ## Get clusters in order
  clusters <- cluster_vec %>% unique()
  cluster_vec <- factor(as.character(cluster_vec), levels = clusters[order(as.numeric(extractClusterNumber(clusters)))])
  return(cluster_vec)

}

extractClusterNumber <- function(strs){

  p <- NULL
  i <- 1
  for(str1 in strs){
    p[[i]] <- strsplit(str1, "[ ]")[[1]][1]
    i <- i + 1
  }

  return(p)

}

reorderClusters <- function(cluster_vec){

  ## Get clusters in order
  clusters <- cluster_vec %>% unique()
  cluster_vec <- factor(as.character(cluster_vec), levels = clusters[order(as.numeric(extractClusterNumber(clusters)))])
  return(cluster_vec)

}

facets_nice <- theme(strip.background = element_rect(fill="grey96"), strip.text = element_text(colour = 'black'))
add_guide   <- guides(colour = guide_legend(override.aes = list(size=5)))

getPalette  <- colorRampPalette(brewer.pal(9, "Set3"))
getPalette2 <- colorRampPalette(brewer.pal(8, "Set2"))
getPalette3 <- colorRampPalette(brewer.pal(9, "Set1"))

extractClusterNumber <- function(strs){

  p <- NULL
  i <- 1
  for(str1 in strs){
    p[[i]] <- strsplit(str1, "[ ]")[[1]][1]
    i <- i + 1
  }

  return(p)

}

getLatentUMAP <- function(seurat_object){

  umap_df           <- seurat_object[["latent"]]@cell.embeddings %>% uwot::umap()
  colnames(umap_df) <- c("latent_umap1", "latent_umap2")
  rownames(umap_df) <- colnames(seurat_object)
  umap_df           <- CreateDimReducObject(key = "latent_umap", embeddings = as.matrix(x = umap_df))
  seurat_object[['latent_umap']] <- umap_df

  return(seurat_object)

}

plotLatentUmap <- function(viz_df, cluster){

  viz_df_temp <- data.frame(viz_df, "seurat_cluster" = cluster)
  nClusters   <- unique(cluster) %>% length

  ## Visualise
  umap_mean <- data.frame(aggregate(umap_1 ~ seurat_cluster, viz_df_temp, median), umap_2 = aggregate(umap_2 ~ seurat_cluster, viz_df_temp, median)[,2])


  ## Plot UMAPs with TCRGP predictions highlighted
  ggplot() +
    geom_point(data = viz_df_temp, aes(x = umap_1, y = umap_2, color = seurat_cluster), size = 0.8) +

    # stat_ellipse(data = viz_df_temp, geom = "polygon", aes(x = umap_1, y = umap_2, color = seurat_cluster, fill = seurat_cluster), alpha = 0.1, lty = "dotted") +
    ggrepel::geom_label_repel(data = umap_mean, aes(x = umap_1, y = umap_2, color = seurat_cluster, label = seurat_cluster), size = 5, color = "black") +

    theme_void() + theme(legend.position = "none") +
    scale_color_manual(values = getPalette(nClusters)) +
    scale_fill_manual(values = getPalette(nClusters)) + labs()


}

extractName = function(str1){
  # strsplit(str1, "[_]")[[1]][1]
  sub("\\_.*", "", str1)
}

extractFileName = function(str1){
  # strsplit(str1, "[_]")[[1]][1]
  sub(".*\\/", "", str1)
}

extractSeuratName <- function(str1){

  strsplit(str1, "[//]")[[1]][4]

}

plotQcViolin <- function(viz_df, var_to_plot, grouping, min, max){

  ## Plot univariate violin plots with filter thresholds

  # @ params:
  # viz_df = df that contains qc-analysis results and covariates of interest
  # var_to_plot = char, a column name that contains the variable to plot
  # grouping = char, a column name that contains the x-axis grouping
  # min = num, min value for variable
  # max = num, max value for variable

  viz_df_temp <- viz_df %>% dplyr::select(var_to_plot)

  label_df_min <- ifelse(viz_df_temp > min, "above", "below") %>% table
  label_df_max <- ifelse(viz_df_temp < max, "above", "below") %>% table

  ggplot(data = viz_df, aes_string(x = grouping, y = var_to_plot, fill = grouping)) +
    geom_violin(alpha = 0.5) +
    # geom_boxplot(width = 0.1, alpha = 0.8, outlier.shape = NA) +

    geom_hline(yintercept = min, linetype = "dotted") +
    geom_hline(yintercept = max, linetype = "dotted") +

    annotate(geom = "text", x = 2.5, y = min, label = paste("Below the line:\n", label_df_min[2]), fontface = "italic") +
    annotate(geom = "text", x = 2.5, y = max, label = paste("Above the line:\n", label_df_max[2]), fontface = "italic") +

    labs(x = "", title = var_to_plot) + theme(legend.position = "none")

}






# getClusterPhenotypes <- function(clusters){
#
#   clusters <- plyr::revalue(clusters, replace = c(
#
#     "0"  = "00_NK_NK",
#     "1"  = "01_CD4_naive",
#     "2"  = "02_CD4_EM/Th1-like",
#     "3"  = "03_CD8_EM",
#     "4"  = "04_CD8_TRM",
#     "5"  = "05_CD4CD8_naive",
#     "6"  = "06_CD4CD8_naive",
#     "7"  = "07_CD4_TRM/Th2",
#     "8"  = "08_other_MAIT",
#     "9"  = "09_B",
#     "10" = "10_CD8_effector/exhausted",
#     "11" = "11_CD4_Treg",
#     "12" = "12_B_B/T",
#     "13" = "13_B_?",
#     "14" = "14_NK_?",
#     "15" = "15_B_?",
#     "16" = "16_CD8_?",
#     "17" = "17_monocyte_promonocyte",
#     "18" = "18_monocyte_monocyte",
#     "19" = "19_NK_?",
#     "20" = "20_other_HSC-like?",
#     "21" = "21_other_pDC",
#     "22" = "22_other_junk",
#     "23" = "23_other_plasma"))
#
#
#   return(clusters)
#
#
# }



getClusterPhenotypes <- function(clusters){

  clusters <- plyr::revalue(clusters, replace = c(

    "0"  = "0 CD8 effector IFNg",
    "1"  = "1 LGLL FM1656" ,
    "2"  = "2 CD4 EM" ,
    "3"  = "3 Monocytes" ,
    "4"  = "4 CD8 effector" ,
    "5"  = "5 CD4 naive/CM" ,
    "6"  = "6 LGLL low-TCR FM1656 specific" ,
    "7"  = "7 CD8 EM" ,
    "8"  = "8 NK CD56bright" ,
    "9"  = "9 LGLL low-TCR FM2481 specific" ,
    "10" = "10 LGLL low-TCR FM1656 specific" ,
    "11" = "11 LGLL FM1656" ,
    "12" = "12 LGLL FM805" ,
    "13" = "13 Monocytes" ,
    "14" = "14 B-cell" ,
    "15" = "15 CD8 effector" ,
    "16" = "16 CD4 naive/CM" ,
    "17" = "17 Doublets" ,
    "18" = "18 CD4 tregs" ,
    "19" = "19 CD8 effector" ,
    "20" = "20 CD4" ,
    "21" = "21 lymphocytes" ,
    "22" = "22 CD4",
    "23" = "23 cycling",
    "24" = "24 Monocytes",
    "25" = "25 pDC"))

  return(clusters)

}


## from dipa
getClusterPhenotypes <- function(clusters){

  clusters <- plyr::revalue(clusters, replace = c(

    "0"  = "0 LGLL",
    "1"  = "1 LGLL" ,
    "2"  = "2 CD4 EM" ,
    "3"  = "3 Monocyte CD16-" ,
    "4"  = "4 LGLL" ,
    "5"  = "5 CD4 naive/CM" ,
    "6"  = "6 LGLL" ,
    "7"  = "7 CD8 EM" ,
    "8"  = "8 NK CD56bright" ,
    "9"  = "9 CD8 effector" ,
    "10" = "10 CD8 effector" ,
    "11" = "11 LGLL" ,
    "12" = "12 LGLL" ,
    "13" = "13 Monocyte CD16+" ,
    "14" = "14 B-cell" ,
    "15" = "15 LGLL" ,
    "16" = "16 CD8" ,
    "17" = "17 Doublets (T/monocyte)" ,
    "18" = "18 CD4 tregs" ,
    "19" = "19 CD8 effector" ,
    "20" = "20 CD4 emra" ,
    "21" = "21 CD8" ,
    "22" = "22 CD4 activated",
    "23" = "23 cycling",
    "24" = "24 cDC",
    "25" = "25 pDC"))

  return(clusters)

}




getClusterPhenotypesPublic <- function(clusters){

  clusters <- plyr::revalue(clusters, replace = c(

    "0"  = "0 CD8 effector",
    "1"  = "1 CD8 effector IFNg" ,
    "2"  = "2 Monocytes 1" ,
    "3"  = "3 CD4" ,
    "4"  = "4 LGLL FM1656" ,
    "5"  = "5 CLL" ,
    "6"  = "6 CD4 naive/CM" ,
    "7"  = "7 NK CD56bright?" ,
    "8"  = "8 B-cell" ,
    "9"  = "9 LGLL low-TCR FM1656 specific" ,
    "10" = "10 LGLL low-TCR FM1656 specific" ,
    "11" = "11 Monocytes 2" ,
    "12" = "12 CD4 naive/CM" ,
    "13" = "13 Doublets" ,
    "14" = "14 LGLL FM805" ,
    "15" = "15 CD8 effector" ,
    "16" = "16 CD4" ,
    "17" = "17 Monocytes 3" ,
    "18" = "18 Monocytes 4" ,
    "19" = "19 CD8 EM" ,
    "20" = "20 pDC" ,
    "21" = "21 cycling" ,
    "22" = "22 ",
    "23" = "23 "))

  return(clusters)

}


# getClusterPhenotypesPublicImmune <- function(clusters){
#
#   clusters <- plyr::revalue(clusters, replace = c(
#
#     "0"  = "0 CD8 effector",
#     "1"  = "1 Monocytes 1" ,
#     "2"  = "2 CD4" ,
#     "3"  = "3 CD8 effector IFNg" ,
#     "4"  = "4 NK CD56bright?" ,
#     "5"  = "5 CD4 naive/CM" ,
#     "6"  = "6 LGLL FM1656" ,
#     "7"  = "7 Monocytes 2" ,
#     "8"  = "8 B-cell 1" ,
#     "9"  = "9 B-cell 2" ,
#     "10" = "10 CD4 naive/CM" ,
#     "11" = "11 Doublets" ,
#     "12" = "12 CD8 effector" ,
#     "13" = "13 LGLL low-TCR FM1656 specific" ,
#     "14" = "14 CD4" ,
#     "15" = "15 CD4" ,
#     "16" = "16 Monocytes 3" ,
#     "17" = "17 Monocytes 4" ,
#     "18" = "18 pDC" ,
#     "19" = "19 cycling" ,
#     "20" = "20 " ))
#
#   return(clusters)
#
# }

# getClusterPhenotypesPublicImmune <- function(clusters){
#
#   clusters <- plyr::revalue(clusters, replace = c(
#
#     "0"  = "0 Myeloid",
#     "1"  = "1 CD8 effector" ,
#     "2"  = "2 CD4 Treg" ,
#     "3"  = "3 CD8 EM IFNg" ,
#     "4"  = "4 NK" ,
#     "5"  = "5 CD8 effector/exhausted" ,
#     "6"  = "6 CD4 naive/CM" ,
#     "7"  = "7 FM1656 specific" ,
#     "8"  = "8 Myeloid" ,
#     "9"  = "9 FM1656 specific" ,
#     "10" = "10 low quality" ,
#     "11" = "11 B-cell" ,
#     "12" = "12 B-cell" ,
#     "13" = "13 " ,
#     "14" = "14 CD8" ,
#     "15" = "15 Doublets" ,
#     "16" = "16 CD8" ,
#     "17" = "17 CD8" ,
#     "18" = "18 CD8" ,
#     "19" = "19 CD4" ,
#     "20" = "20 FM1656 specific" ,
#     "21" = "21 Myeloid" ,
#     "22" = "22 cDC",
#     "23" = "23 pDC",
#     "24" = "24 cycling",
#     "25" = "25 ",
#     "26" = "26 ",
#     "27" = "27 "
#
#
#     ))
#
#   return(clusters)
#
# }

getClusterPhenotypesPublicImmune <- function(clusters){

  clusters <- plyr::revalue(clusters, replace = c(

    "0"  = "0 CD8 IFNg",
    "1"  = "1 CD4" ,
    "2"  = "2 Myeloid" ,
    "3"  = "3 NK CD56dim" ,
    "4"  = "4 B-cells" ,
    "5"  = "5 CD8 effector" ,
    "6"  = "6 CD8 effector" ,
    "7"  = "7 CD8 effector" ,
    "8"  = "8 Myeloid" ,
    "9"  = "9 CD8 naive/CM" ,
    "10" = "10 CD4/myeloid doublets" ,
    "11" = "11 CD8" ,
    "12" = "12 gd" ,
    "13" = "13 CD8" ,
    "14" = "14 unspecific" ,
    "15" = "15 cDC" ,
    "16" = "16 cycling" ,
    "17" = "17 pDC" ,
    "18" = "18 Myeloid" ,
    "19" = "19 B-cells plasma cells"))

  return(clusters)

}




getClusterPhenotypesPublicImmuneNonLGLL <- function(clusters){

  clusters <- plyr::revalue(clusters, replace = c(

    "0"  = "0 CD8 EM IFNg LAG3",
    "1"  = "1 CD4 naive/CM" ,
    "2"  = "2 Monocytes CD16-" ,
    "3"  = "3 CD8" ,
    "4"  = "4 NK" ,
    "5"  = "5 B-cell" ,
    "6"  = "6 Monocytes CD16+" ,
    "7"  = "7 low quality T cells" ,
    "8"  = "8 CD4 naive/CM" ,
    "9"  = "9 Monocytes CD16-" ,
    "10" = "10 CD8" ,
    "11" = "11 gd" ,
    "12" = "12 CD8" ,
    "13" = "13 CD4" ,
    "14" = "14 cDC" ,
    "15" = "15 cycling" ,
    "16" = "16 pDC" ,
    "17" = "17 doublets" ,
    "18" = "18 " ,
    "19" = "19 B-cell plasma" ))

  return(clusters)

}




getClusterPhenotypesLGL <- function(clusters){

  clusters <- plyr::revalue(clusters, replace = c(

    "0"  = "0 GZMB+ STAT3+ TIGIT+ LAG3+",
    "1"  = "1 GAPDH+ MT2A+" ,
    "2"  = "2 CD3 low CD8 low TOB1+" ,
    "3"  = "3 CD3 low CD8 low IL2RB+" ,
    "4"  = "4 highly cytotoxic KLRB1+ FCGR3A+" ,
    "5"  = "5 EM GZMK+ STAT3+" ,
    "6"  = "6 PRF1+ KLRG1+" ,
    "7"  = "7 CD74+" ,
    "8"  = "8 IFNG+" ,
    "9"  = "9 most cytotoxic PDCD1+ CXCR3+" ,
    "10" = "10 RUNX3+" ,
    "11" = "11 " ,
    "12" = "12 STAT3+" ,
    "13" = "13 naive"))

  return(clusters)

}


getClusterPhenotypesLGL <- function(clusters){

  plyr::revalue(clusters, replace = c( "0"  = "0 PRF11+ LAG3+ IFNG+",
                                       "1"  = "1 GZMA+ KLRB1+" ,
                                       "2"  = "2 GNLY+ NKG7+" ,
                                       "3"  = "3 STAT3wt IFNG+" ,
                                       "4"  = "4 STAT3wt FCGR3A + CXCR3+" ,
                                       "5"  = "5 GZMK+" ,
                                       "6"  = "6 KLRB1+" ,
                                       "7"  = "7 GZMK+ CXCR3+" ,
                                       "8"  = "8 IFNG+ TIGIT+" ,
                                       "9"  = "9 STAT3wt PDCD1+" ,
                                       "10" = "10 " ,
                                       "11" = "11 CCL5+" ,
                                       "12" = "12 STAT3wt CD4+" ,
                                       "13" = "13 STAT3wt PRF1+" ))
}


getClusterPhenotypesLGL <- function(clusters){

  plyr::revalue(clusters, replace = c( "0"  = "0 PRF11+ LAG3+ IFNG+",
                                       "1"  = "1 GZMA+ KLRB1+" ,
                                       "2"  = "2 GNLY+ NKG7+" ,
                                       "3"  = "3 STAT3wt IFNG+" ,
                                       "4"  = "4 STAT3wt FCGR3A + CXCR3+" ,
                                       "5"  = "5 GZMK+" ,
                                       "6"  = "6 KLRB1+" ,
                                       "7"  = "7 GZMK+ CXCR3+" ,
                                       "8"  = "8 IFNG+ TIGIT+" ,
                                       "9"  = "9 STAT3wt PDCD1+" ,
                                       "10" = "10 " ,
                                       "11" = "11 CCL5+" ,
                                       "12" = "12 STAT3wt CD4+" ,
                                       "13" = "13 STAT3wt PRF1+" ))
}



#
# getClusterPhenotypesLGL <- function(clusters){
#
#   plyr::revalue(clusters, replace = c("0"	 = "0 PRF1+ LAG3+ CCL3+ STAT3+",
#                                       "1"	 = "1 FM1656 1",
#                                       "2"	 = "2 FM2481 IL2RB+ NFKBIA+",
#                                       "3"	 = "3 FM805 GZMK+ STAT3+",
#                                       "4"	 = "4 IFNg+ TIGIT+",
#                                       "5"	 = "5 ",
#                                       "6"	 = "6 Uncertain LGLL phenotype",
#                                       "7"	 = "7 CTLA4+ STAT3+"))
# }

getClusterPhenotypesXXX <- function(clusters){

  clusters <- plyr::revalue(clusters, replace = c(

    "0"  = "0 ",
    "1"  = "1 " ,
    "2"  = "2 " ,
    "3"  = "3 " ,
    "4"  = "4 " ,
    "5"  = "5 " ,
    "6"  = "6 " ,
    "7"  = "7 " ,
    "8"  = "8 " ,
    "9"  = "9 " ,
    "10" = "10 " ,
    "11" = "11 " ,
    "12" = "12 " ,
    "13" = "13 " ,
    "14" = "14 " ,
    "15" = "15 " ,
    "16" = "16 " ,
    "17" = "17 " ,
    "18" = "18 " ,
    "19" = "19 " ,
    "20" = "20 " ,
    "21" = "21 " ,
    "22" = "22 ",
    "23" = "23 ",
    "24" = "24 ",
    "25" = "25 "))

  return(clusters)

}







getClusterCoarsePhenotype <- function(clusters){

  # func_temp <- function(str1){strsplit(clusters, "[ ]")[[1]][2]}
  # vec <- lapply(clusters, func_temp)
  # vec <- do.call(what = "c", vec)
  # return(vec)

  clusters <- substr(clusters, 1, 2) %>% as.numeric() %>% as.character()

  clusters <- plyr::revalue(clusters, replace = c(

    "0"  = "NK",
    "1"  = "CD4",
    "2"  = "CD4",
    "3"  = "CD8",
    "4"  = "CD8",
    "5"  = "CD4CD8",
    "6"  = "CD4CD8",
    "7"  = "CD4",
    "8"  = "other",
    "9"  = "B",
    "10" = "CD8",
    "11" = "CD4",
    "12" = "B",
    "13" = "B",
    "14" = "NK",
    "15" = "B",
    "16" = "CD4",
    "17" = "other",
    "18" = "other",
    "19" = "NK",
    "20" = "CD4CD8",
    "21" = "other",
    "22" = "other",
    "23" = "other"))

  return(clusters)

}





## Helper function to break TCRb names into clinical data
breakName <- function(name_temp){

  ## Breaks the filename into relevant information

  # "2121_PR_PB_LAG3_R_0m_MNC.txt"

  name      <- substr(name_temp, 1, 4)
  response  <- substr(name_temp, 6, 7)
  type      <- substr(name_temp, 9, 10)
  regimen   <- substr(name_temp, 12, 15)
  overall   <- substr(name_temp, 17, 17)
  timepoint <- substr(name_temp, 19, 20)
  celltype  <- substr(name_temp, 22, 24)
  io_stat   <- substr(name_temp, 26, 27)


  for(i in 1:length(type)){
    if(type[i] == "PB"){type[i] = "Blood"}
    if(type[i] == "TX"){type[i] = "Tumor"}
  }

  for(i in 1:length(regimen)){
    if(regimen[i] == "CTLA"){regimen[i] = "antiCTLA4"}
    if(regimen[i] == "iPD1"){regimen[i] = "antiPD1"}
    if(regimen[i] == "NiIp"){regimen[i] = "antiPD1+antiCTLA4"}
    if(regimen[i] == "IpNi"){regimen[i] = "antiCTLA4+antiPD1"}
    if(regimen[i] == "LAG3"){regimen[i] = "antiPD1+antiLAG3"}
    if(regimen[i] == "NANA"){regimen[i] = "TreatmentNA"}
  }

  for(i in 1:length(overall)){
    if(overall[i] == "U"){overall[i] = "NA"}
  }

  for(i in 1:length(io_stat)){
    if(io_stat[i] == "pr"){io_stat[i] = "Prior.IO"}
    if(io_stat[i] == "nv"){io_stat[i] = "IO.naive"}

  }

  return(data.frame(name, response, type, regimen, overall, timepoint, celltype, io_stat, stringsAsFactors = F))

}




fixSeurat <- function(seurat_object){

  ## Fix meta data if it brokes

  meta.data           <- seurat_object@meta.data
  count.data          <- seurat_object@assays$RNA@counts
  scale.data          <- seurat_object@assays$RNA@scale.data
  # hvg                 <- VariableFeatures(seurat_object)

  # pca_dimred          <- seurat_object[["pca"]]
  # umap_dimred         <- seurat_object[["umap"]]
  latent_dimred       <- seurat_object[["latent"]]
  latent_umap_dimred  <- seurat_object[["latent_umap"]]

  rownames(meta.data) <- meta.data$barcode

  old_idents <- Idents(seurat_object)
  new_seurat <- CreateSeuratObject(counts = count.data)

  new_seurat@meta.data             <- meta.data
  new_seurat@assays$RNA@counts     <- count.data
  new_seurat@assays$RNA@scale.data <- scale.data
  # VariableFeatures(seurat_object)  <- hvg

  # new_seurat[["pca"]]              <- pca_dimred
  # new_seurat[["umap"]]             <- umap_dimred
  new_seurat[["latent"]]           <- latent_dimred
  new_seurat[["latent_umap"]]      <- latent_umap_dimred
  Idents(new_seurat) <- old_idents
  return(new_seurat)

}





extractCoarsePhenotype <- function(strs){

  p <- NULL
  i <- 1
  for(str1 in strs){
    p[[i]] <- strsplit(str1, "[ ]")[[1]][2]
    i <- i + 1
  }

  return(p)

}
