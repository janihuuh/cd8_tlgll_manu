
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

