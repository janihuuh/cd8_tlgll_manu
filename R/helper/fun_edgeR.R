
build_edgeR_objects <- function(counts, group1, group2, group1_name = "A", group2_name = "B"){
  
  require(edgeR)
  
  # in: counts-df, vector of samples belonging to group 1 and 2;
  # names for groups 1 and 2
  # out: edgeR-object, counted with common and tagwise dispersions
  
  
  ## Add design, e.g. which group compare to which

  keep <- c(group1, group2)
  design <- rep(c(group1_name, group2_name), c(length(group1), length(group2)))

  ## Make edgeR-compatible object
  cds_temp <- DGEList(counts = counts[ ,keep], group = design)
    
  # Remove lowly expressed genees
  cds_temp <- cds_temp[rowSums(1e+06 * cds_temp$counts/expandAsMatrix(cds_temp$samples$lib.size, dim(cds_temp)) > 1) >= 3, ]
  
  ## Add prior to tagwise-dispersion, i.e. how much dispersion is 
  # "smoothened" towards common dispersion (see vignette for details)
  # Prior: by recommendation or by heuristic
  # prior <- floor(50 / length(keep) - 2)
  prior = 10
  
  # Count normalisation factors and common and tagwise dispersions
  cds_temp <- calcNormFactors(cds_temp)
  cds_temp <- estimateCommonDisp(cds_temp)
  cds_temp <- estimateTagwiseDisp(cds_temp, prior.n = prior)
  
  return(cds_temp)
  
}



get_bm <- function(genes){
  
  require(biomaRt)
  ## Anotate ensembl-genes
  
  # in: vector of ensembl_genes
  # out: df with annotation for ensembl genes
  ensembl    = useMart("ENSEMBL_MART_ENSEMBL")
  ensembl    = useDataset("hsapiens_gene_ensembl", mart = ensembl)
  attributes = c("ensembl_gene_id", "external_gene_name", "description", "gene_biotype")
  
  gene.df <- getBM(attributes = attributes,
                   filters = "ensembl_gene_id",
                   values = genes,
                   mart = ensembl)
  
}




estimate_de <- function(cds_temp, group1_name = "A", group2_name = "B"){
  
  
  # in: edgerR object, names for groups 1 and 2
  # out: DE-gene df, with annotations from biomaRt
  

  ## Do the DE-test
  de.cmn <- exactTest(cds_temp, dispersion = "common",  pair = c(group2_name, group1_name))$table
  de.tgw <- exactTest(cds_temp, dispersion = "tagwise", pair = c(group2_name ,group1_name))$table
  
  ## Correct p-values with Benjamini-Hochberg
  de.cmn$qval <- p.adjust(de.cmn$PValue, method = "BH")
  de.tgw$qval <- p.adjust(de.tgw$PValue, method = "BH")
  
  ## Annotate significance and direction of genes
  de.cmn$sigf <- ifelse(de.cmn$qval < 0.05 & abs(de.cmn$logFC) > 1, "Significant", "Not_significant")
  de.tgw$sigf <- ifelse(de.tgw$qval < 0.05 & abs(de.tgw$logFC) > 1, "Significant", "Not_significant")

  de.cmn$direction <- ifelse(de.cmn$logFC > 0, "Up", "Down")
  de.tgw$direction <- ifelse(de.tgw$logFC > 0, "Up", "Down")
  
  print("De.cmn")
  print(table(de.cmn$sigf))

  print("De.tgq")
  print(table(de.tgw$sigf))

  ## Merge different dispersions
  colnames(de.cmn) <- paste0(colnames(de.cmn), "_cmn")
  colnames(de.tgw) <- paste0(colnames(de.tgw), "_tgw")
  
  de_df <- cbind(de.cmn, de.tgw)
  
  
  ## Add biomart annotation
  de_bm <- get_bm(rownames(de_df))
  de_df_full <- merge(de_df, de_bm, by.x = 0, by.y = "ensembl_gene_id")
  
  print(paste("Lost information of", nrow(de_df_full) - nrow(de_df), "genes"))
  
  return(de_df_full)
  
}











enrich_hypergeometric <- function(de_df_full, direction, universe_df){
  
  require(clusterProfiler)
  
  # Enrichment in hypergeometric test
  # in: de_df with BM-annotations, direction (Up/Down) and universe to count the enrichment
  # from, e.g. all genes available (~60 000) or all genes expressed in the data (usually ~20 000)
    
  # out: df with enrichment results
  
  df_cmn      <- subset(de_df_full, sigf_cmn == "Significant" & direction_cmn == direction)
  sigf_cmn    <- df_cmn[order(df_cmn$qval_cmn, decreasing = T), ]$external_gene_name
  enrich_cmn  <- enricher(sigf_cmn, universe = universe_df$external_gene_name, TERM2GENE = hallmark)
  results_cmn <- do.call(rbind, enrich_cmn@result)
  results_cmn <- as.data.frame(t(results_cmn))

  df_tgw      <- subset(de_df_full, sigf_tgw == "Significant" & direction_tgw == direction)
  sigf_tgw    <- df_tgw[order(df_tgw$qval_tgw, decreasing = T), ]$external_gene_name
  enrich_tgw  <- enricher(sigf_tgw, universe = universe_df$external_gene_name, TERM2GENE = hallmark)
  results_tgw <- do.call(rbind, enrich_tgw@result)
  results_tgw <- as.data.frame(t(results_tgw))

  results_cmn$dispersion_type <- "Common"
  results_tgw$dispersion_type <- "Tagwise"
  results                     <- rbind(results_cmn, results_tgw)
  
  colnames(results)
  results[,c(5:7, 9)] <- sapply(results[,c(5:7, 9)], function(x) {as.numeric(as.character(x))})
  return(results)
  
}

