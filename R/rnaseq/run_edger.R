
# Datamanipulation and plotting
library(gridExtra)
library(ggplot2)
library(reshape2)
library(grid)


# Analyses
library(edgeR)
library(biomaRt)
library(org.Hs.eg.db)
library(clusterProfiler)
library(DESeq)

theme_set(theme_classic())


## Slight preprocssing of data
lgl_counts     <- fread("data/RNAseq/lgl_counts_new.txt")
genenames      <- lgl_counts$gene
lgl_counts     <- lgl_counts[,-c(1:2)] %>% as.data.frame()
rownames(lgl_counts) <- genenames

patients       <- fread("data/clinical/lgl_patients_corrected_full_jani.txt")
patients$STAT3 <- patients$Stat3_status
patients$ID    <- patients$FM

## Add different study designs

# 1) STAT3mt all vs healthy
# 2) STAT3mt all vs STAT3wt LGL

# 3) STAT3mt CD8 and CD4/CD8 vs healthy CD8
# 4) STAT3mt CD8 and CD4/CD8 vs STAT3wt LGL

# 5) STAT3mt CD8 vs healthy CD8
# 6) STAT3mt CD8 vs STAT3wt LGL

# and bonus;
# 7) Healthy CD8 vs healthy CD4

stat3mt_all    <- which(patients$STAT3 != "WT")
stat3mt_cd8    <- which(patients$STAT3 != "WT" & patients$Cell_type == "CD8")
stat3mt_cd4    <- which(patients$STAT3 != "WT" & patients$Cell_type == "CD4")
stat3mt_cd4cd8 <- which(patients$STAT3 != "WT" & patients$Cell_type == "CD4CD8")

healthy_all    <- which(patients$Status == "Healthy")
healthy_cd8    <- which(patients$Status == "Healthy" & patients$Cell_type == "CD8")
healthy_cd4    <- which(patients$Status == "Healthy" & patients$Cell_type == "CD4")

lgl_STAT3_wt  <-  which(patients$STAT3 == "WT" & patients$Status == "T-LGL" & patients$Cell_type == "CD8")


# Graphic titles
a_title = "STAT3MT all vs healthy all"
b_title = "STAT3MT all vs STAT3wt LGL"

c_title = "STAT3MT CD8 and CD4/CD8 vs healthy CD8"
d_title = "STAT3MT CD8 and CD4/CD8 vs STAT3wt LGL"

e_title = "STAT3MT CD8 vs healthy CD8"
f_title = "STAT3MT CD8 vs STAT3wt LGL"

g_title = "Healthy CD8 vs healthy CD4"


a_name <- paste(c(strsplit(a_title, " ")[[1]]), collapse = "_")
b_name <- paste(c(strsplit(b_title, " ")[[1]]), collapse = "_")
c_name <- paste(c(strsplit(c_title, " ")[[1]]), collapse = "_")
d_name <- paste(c(strsplit(d_title, " ")[[1]]), collapse = "_")
e_name <- paste(c(strsplit(e_title, " ")[[1]]), collapse = "_")
f_name <- paste(c(strsplit(f_title, " ")[[1]]), collapse = "_")
g_name <- paste(c(strsplit(g_title, " ")[[1]]), collapse = "_")

c_name <- paste0(substr(c_name, 1, 19), substr(c_name, 21, nchar(c_name)))
d_name <- paste0(substr(c_name, 1, 19), substr(c_name, 21, nchar(c_name)))

name_list <- list(a_name, b_name, c_name, d_name, e_name, f_name, g_name)

grob_to_pdf <- function(grob_temp, title, name = "", folder = "", width = 9, height = 9){

  ggsave(grob_temp, file = paste0(folder, title, "_", name, ".pdf"),  width = width, height = height)

}


## Build diagnostic PCA:s
calc_and_plot_pca <- function(group1, group2, title = ""){

  # group1 = stat3mt_all
  # group2 = healthy_all

  pca_temp     <- prcomp(t(lgl_counts[,c(group1, group2)]))
  patient_temp <- patients[c(group1, group2),]

  a <- ggplot(as.data.frame(pca_temp$x), aes(x = PC1, y = PC2, label = patient_temp$ID, color = patient_temp$STAT3)) +
    geom_point(size = 3) + labs(color = "STAT3", title = "STAT3") + geom_text(nudge_y = 100, size = 3)

  b <- ggplot(as.data.frame(pca_temp$x), aes(x = PC1, y = PC2, label = patient_temp$ID, color = patient_temp$Celltype)) +
    geom_point(size = 3) + labs(color = "Celltype", title = "Celltype") + geom_text(nudge_y = 101, size = 3)

  return(grid.arrange(a,b,top = textGrob(title,gp=gpar(fontsize=20,font=3)), ncol = 2))

}


a <- calc_and_plot_pca(group1 = stat3mt_all, group2 = healthy_all, title = a_title)
b <- calc_and_plot_pca(stat3mt_all, lgl_STAT3_wt, title = b_title)

c <- calc_and_plot_pca(c(stat3mt_cd8, stat3mt_cd4cd8), healthy_cd8, title = c_title)
d <- calc_and_plot_pca(c(stat3mt_cd8, stat3mt_cd4cd8), lgl_STAT3_wt, title = d_title)

e <- calc_and_plot_pca(stat3mt_cd8, healthy_cd8, title = e_title)
f <- calc_and_plot_pca(stat3mt_cd8, lgl_STAT3_wt, title = f_title)

g <- calc_and_plot_pca(healthy_cd8, healthy_cd4, title = g_title)

dir.create("results/PCA/", showWarnings = F)

pdf("results/PCA/pca_meta.pdf", width = 16, height = 38)
grid.arrange(a,b,c,d,e,f,g,ncol = 1)
dev.off()

# Save individual plots
grob_list <- list(a,b,c,d,e,f,g)
for(i in 1:7){
  grob_to_pdf(grob_list[[i]], title = "PCA", name = name_list[[i]], folder = "results/PCA/", width = 10, height = 4.5)
}



# ========== edgeR ===========


build_edgeR_objects <- function(group1, group2, group1_name = "A", group2_name = "B"){

  keep <- c(group1, group2)

  design <- rep(c(group1_name, group2_name), c(length(group1), length(group2)))
  
  cds_temp <- DGEList(counts = lgl_counts[ ,keep], group = design)

  dim(lgl_counts[ ,keep]) == length(design)

  # Remove lowly expressed genees
  cds_temp <- cds_temp[rowSums(1e+06 * cds_temp$counts/expandAsMatrix(cds_temp$samples$lib.size, dim(cds_temp)) > 1) >= 3, ]

  # Prior: by recommendation or by heuristic
  # prior <- floor(50 / length(keep) - 2)
  prior = 10

  # Count normalisation factors and cmn and tgw dispersions
  cds_temp <- calcNormFactors(cds_temp)
  cds_temp <- estimateCommonDisp(cds_temp)
  cds_temp <- estimateTagwiseDisp(cds_temp, prior.n = prior)

  return(cds_temp)

}

plot_mean_var <- function(cds_temp, title = ""){

  plotMeanVar(cds_temp ,show.raw.vars=TRUE , show.tagwise.vars=TRUE , show.binned.common.disp.vars=FALSE , show.ave.raw.vars=FALSE , dispersion.method = "qcml" , NBline = TRUE , nbins = 100 , #these are arguments about what is plotted
              pch = 16 , xlab ="Mean Expression (Log10 Scale)" , ylab = "Variance (Log10 Scale)" , main = paste0("Mean-Variance Plot \n", title))


}

estimate_de <- function(cds_temp, group1_name = "A", group2_name = "B"){

  # Estimate de-genes
  # group1_name = "STAT3mt_all"
  # group2_name = "Healthy_all"
  # cds_temp = stat3mt_all_healthy_all_cds

  de.cmn <- exactTest(cds_temp, dispersion = "common",  pair = c(group2_name, group1_name))$table
  de.tgw <- exactTest(cds_temp, dispersion = "tagwise", pair = c(group2_name ,group1_name))$table

  de.cmn$qval <- p.adjust(de.cmn$PValue, method = "BH")
  de.tgw$qval <- p.adjust(de.tgw$PValue, method = "BH")

  de.cmn$sigf <- ifelse(de.cmn$qval < 0.05 & abs(de.cmn$logFC) > 1, "Significant", "Not_significant")
  de.tgw$sigf <- ifelse(de.tgw$qval < 0.05 & abs(de.tgw$logFC) > 1, "Significant", "Not_significant")

  de.cmn$direction <- ifelse(de.cmn$logFC > 0, "Up", "Down")
  de.tgw$direction <- ifelse(de.tgw$logFC > 0, "Up", "Down")

  print("De.cmn")
  print(table(de.cmn$sigf))

  print("De.tgq")
  print(table(de.tgw$sigf))

  colnames(de.cmn) <- paste0(colnames(de.cmn), "_cmn")
  colnames(de.tgw) <- paste0(colnames(de.tgw), "_tgw")

  de_df <- cbind(de.cmn, de.tgw)

  xprs  <- data.frame(Mean_exprs = rowSums(cds_temp$counts))
  de_df <- merge(de_df, xprs, by = 0)
  rownames(cds_temp$counts)
  return(de_df)
  
  ## Add biomart annotation
  # de_bm <- get_bm(rownames(de_df))

#  de_df_full <- merge(de_df, de_bm, by.x = "Row.names", by.y = "ensembl_gene_id")
  # de_df_full <- merge(de_df, de_bm, by.x = 0, by.y = "ensembl_gene_id")

  # print(paste("Lost information of", nrow(de_df_full) - nrow(de_df), "genes"))

  # return(de_df)

}


get_bm <- function(genes){

  ensembl    = useMart("ENSEMBL_MART_ENSEMBL")
  ensembl    = useDataset("hsapiens_gene_ensembl", mart = ensembl)
  attributes = c("ensembl_gene_id", "external_gene_name", "description", "gene_biotype")

  gene.df <- getBM(attributes = attributes,
                   filters = "ensembl_gene_id",
                   values = genes,
                   mart = ensembl)

}



plot_volcano_plots <- function(de_df, title = ""){

  a <- ggplot(de_df, aes(x = logFC_cmn, y = -log10(PValue_cmn), color = sigf_cmn, size = logCPM_cmn)) + geom_point(alpha = 0.5) + labs(color = "qval < 0.05 and |logFC| > 1", subtitle = "Common") + scale_color_manual(values = c("grey", "dodgerblue")) + theme(legend.position = "top") +
    ggrepel::geom_text_repel(data = subset(de_df, sigf_cmn == "Significant"), aes(x = logFC_cmn, y = -log10(PValue_cmn), label = Row.names), nudge_y = 0.5, fontface = "italic", force = 3)

  b <- ggplot(de_df, aes(x = logFC_tgw, y = -log10(PValue_tgw), color = sigf_tgw, size = logCPM_tgw)) + geom_point(alpha = 0.5) + labs(color = "qval < 0.05 and |logFC| > 1", subtitle = "Tagwise") + scale_color_manual(values = c("grey", "salmon")) + theme(legend.position = "top") +
    ggrepel::geom_text_repel(data = subset(de_df, sigf_tgw == "Significant"), aes(x = logFC_tgw, y = -log10(PValue_tgw), label = Row.names), nudge_y = 0.5, fontface = "italic", force = 3)

  return(grid.arrange(a,b,ncol=2, top = textGrob(title,gp=gpar(fontsize=20,font=3))))

}



# Building edgeR-objects
stat3mt_all_healthy_all_cds      <- build_edgeR_objects(stat3mt_all, healthy_all,   "STAT3mt_all", "Healthy_all")
stat3mt_all_lgl_wtl_cds          <- build_edgeR_objects(stat3mt_all, lgl_STAT3_wt, "STAT3mt_all", "STAT3wt")

stat3mt_cd8_cd48_healthy_all_cds <- build_edgeR_objects(c(stat3mt_cd8, stat3mt_cd4cd8), healthy_cd8, group1_name = "STAT3mt_cd8_cd48", group2_name = "Healthy_cd8")
stat3mt_cd8_cd48_lgl_wt_cds      <- build_edgeR_objects(c(stat3mt_cd8, lgl_STAT3_wt), healthy_cd8, "STAT3mt_cd8_cd48", "STAT3wt")

stat3mt_cd8_healthy_cd8_cds      <- build_edgeR_objects(stat3mt_cd8, healthy_cd8, "STAT3mt_cd8", "Healthy_cd8")
stat3mt_cd8_lgl_wt_cds           <- build_edgeR_objects(stat3mt_cd8, lgl_STAT3_wt, "STAT3mt_cd8", "STAT3wt")

healthy_cd8_healthy_cd4_cds      <- build_edgeR_objects(healthy_cd8, healthy_cd4, "Healthy_cd8", "Healthy_cd4")



# Plot mean-var plots
pdf("results/DE_tables/edgeR/mean_var_plots.pdf", width = 12, height = 20)
par(mfrow=c(5,2))
plot_mean_var(stat3mt_all_healthy_all_cds, a_title)
plot_mean_var(stat3mt_all_lgl_wtl_cds, b_title)

plot_mean_var(stat3mt_cd8_cd48_healthy_all_cds, c_title)
plot_mean_var(stat3mt_cd8_cd48_lgl_wt_cds, d_title)

plot_mean_var(stat3mt_cd8_healthy_cd8_cds, e_title)
plot_mean_var(stat3mt_cd8_lgl_wt_cds, f_title)

plot_mean_var(healthy_cd8_healthy_cd4_cds, g_title)

par(mfrow=c(1,1))
dev.off()



# Estimate DE
stat3mt_all_healthy_all_df      <- estimate_de(cds_temp = stat3mt_all_healthy_all_cds, group1_name = "STAT3mt_all", group2_name = "Healthy_all")
stat3mt_all_lgl_wtl_df          <- estimate_de(stat3mt_all_lgl_wtl_cds, "STAT3mt_all", "STAT3wt")

stat3mt_cd8_cd48_healthy_all_df <- estimate_de(stat3mt_cd8_cd48_healthy_all_cds, "STAT3mt_cd8_cd48", "Healthy_cd8")
stat3mt_cd8_cd48_lgl_wt_df      <- estimate_de(stat3mt_cd8_cd48_lgl_wt_cds, "STAT3mt_cd8_cd48", "STAT3wt")

stat3mt_cd8_healthy_cd8_df      <- estimate_de(stat3mt_cd8_healthy_cd8_cds, "STAT3mt_cd8", "Healthy_cd8")
stat3mt_cd8_lgl_wt_df           <- estimate_de(stat3mt_cd8_lgl_wt_cds, "STAT3mt_cd8", "STAT3wt")

healthy_cd8_healthy_cd4_df      <- estimate_de(healthy_cd8_healthy_cd4_cds, "Healthy_cd8", "Healthy_cd4")


write.table(stat3mt_all_healthy_all_df, paste0("results/DE_tables/edgeR/DE_", a_name, ".txt"), sep = "\t")
write.table(stat3mt_all_lgl_wtl_df, paste0("results/DE_tables/edgeR/DE_", b_name, ".txt"), sep = "\t")

write.table(stat3mt_cd8_cd48_healthy_all_df, paste0("results/DE_tables/edgeR/DE_", c_name, ".txt"), sep = "\t")
write.table(stat3mt_cd8_cd48_lgl_wt_df, paste0("results/DE_tables/edgeR/DE_", d_name, ".txt"), sep = "\t")

write.table(stat3mt_cd8_healthy_cd8_df, paste0("results/DE_tables/edgeR/DE_", e_name, ".txt"), sep = "\t")
write.table(stat3mt_cd8_lgl_wt_df, paste0("results/DE_tables/edgeR/DE_", f_name, ".txt"), sep = "\t")

write.table(healthy_cd8_healthy_cd4_df, paste0("results/DE_tables/edgeR/DE_", g_name, ".txt"), sep = "\t")




# Plot volcanoplots of DE-genes
a <- plot_volcano_plots(stat3mt_all_healthy_all_df, title = a_title)
b <- plot_volcano_plots(stat3mt_all_lgl_wtl_df, title = b_title)

c <- plot_volcano_plots(stat3mt_cd8_cd48_healthy_all_df, title = c_title)
d <- plot_volcano_plots(stat3mt_cd8_cd48_lgl_wt_df,  title = d_title)

e <- plot_volcano_plots(stat3mt_cd8_healthy_cd8_df, title = e_title)
f <- plot_volcano_plots(stat3mt_cd8_lgl_wt_df, title = f_title)

g <- plot_volcano_plots(healthy_cd8_healthy_cd4_df, title = g_title)

pdf("results/DE_tables/edgeR/volcanoplots_meta.pdf", width = 18, height = 38)
grid.arrange(a,b,c,d,e,f,g,ncol=1)
dev.off()



# Save individual plots
grob_list <- list(a,b,c,d,e,f,g)
for(i in 1:7){
  grob_to_pdf(grob_list[[i]], title = "Volcanoplot", name = name_list[[i]], folder = "results/DE_tables/edgeR/volcano/", width = 18, height = 12)
}




# ====== Enrichment analysis

library(GSEABase)
library(grid)


enrich_gsea <- function(df_temp){

  percentage <- length(intersect(df_temp$external_gene_name, hallmark$gene)) / nrow(hallmark)
  print(paste(round(percentage, 3), "% of HALLMARK genes found in df..."))

  # GSEA
  df_cmn          <- df_temp[order(df_temp$PValue_cmn, decreasing = T), ]
  cmn_gsea        <- df_cmn$PValue_cmn
  names(cmn_gsea) <- df_cmn$Row.names
  cmn_gsea        <- cmn_gsea[!duplicated(names(cmn_gsea))]

  gsea_cmn        <- GSEA(cmn_gsea, TERM2GENE = hallmark, verbose = FALSE, pvalueCutoff = 0.05)
  results_cmn     <- do.call(rbind, gsea_cmn@result)
  results_cmn     <- as.data.frame(t(results_cmn))

  df_tgw          <- df_temp[order(df_temp$PValue_tgw, decreasing = T), ]
  tgw_gsea        <- df_tgw$PValue_tgw
  names(tgw_gsea) <- df_tgw$Row.names
  gsea_tgw        <- GSEA(tgw_gsea, TERM2GENE = hallmark, verbose = FALSE)
  results_tgw     <- do.call(rbind, gsea_tgw@result)
  results_tgw     <- as.data.frame(t(results_tgw))

  results <- rbind(results_cmn, results_tgw)

  return(results)


}

enrich_hypergeometric <- function(df_temp, dir_temp){

  # df_temp = stat3mt_cd8_lgl_wt_df
  # dir_temp = "Down"

  # Enrichment as in hypergeometric test
  df_cmn      <- subset(df_temp, sigf_cmn == "Significant" & direction_cmn == dir_temp)
  sigf_cmn    <- df_cmn[order(df_cmn$qval_cmn, decreasing = T), ]$Row.names
  enrich_cmn  <- enricher(sigf_cmn, universe = universe_df, TERM2GENE = hallmark)
  results_cmn <- do.call(rbind, enrich_cmn@result)
  results_cmn <- as.data.frame(t(results_cmn))
  

  df_tgw      <- subset(df_temp, sigf_tgw == "Significant" & direction_tgw == dir_temp)
  sigf_tgw    <- df_tgw[order(df_tgw$qval_tgw, decreasing = T), ]$Row.names
  enrich_tgw  <- enricher(sigf_tgw, universe = universe_df, TERM2GENE = hallmark)
  
  if(!is.null(enrich_tgw)){
    
    results_tgw <- do.call(rbind, enrich_tgw@result)
    results_tgw <- as.data.frame(t(results_tgw))
    
    results_cmn$dispersion_type <- "Common"
    results_tgw$dispersion_type <- "Tagwise"
    results                     <- rbind(results_cmn, results_tgw)
    
    colnames(results)
    results[,c(5:7, 9)] <- sapply(results[,c(5:7, 9)], function(x) {as.numeric(as.character(x))})
    return(results)
    
  }


}

plot_enrichment <- function(enrichment_df, title = ""){

  # enrichment_df = stat3mt_cd8_healthy_cd8_enrich

  enrichment_df$ID <- substr(enrichment_df$ID, 10, nchar(as.character(enrichment_df$ID)))
  interferon_df    <- enrichment_df[grep("*INTERFERON*", as.character(enrichment_df$ID)), ]
  enrichment_df    <- subset(enrichment_df, qvalue < 0.05)

  a <- ggplot(enrichment_df, aes(x = ID, y = Count, fill = dispersion_type)) + geom_bar(stat = "identity", position = "dodge") + coord_flip() + labs(x = "", title = "Significantly (qvalue < 0.05) enriched HALLMARK genesets", fill = "Dispersion") + theme(legend.position = "none")
  b <- ggplot(interferon_df, aes(x = ID, y = qvalue, fill = dispersion_type)) + geom_bar(stat = "identity", position = "dodge") + coord_flip()  + labs(x = "", title = "Interferon HALLMARK genesets", fill = "Dispersion") + geom_hline(yintercept = 0.05, linetype = "dotted")

  return(grid.arrange(a,b,ncol=1, top = textGrob(title,gp=gpar(fontsize=20,font=3))))

}


hallmark   <- read.gmt("/Users/janihuuh/Dropbox/applications/GSEA/h.all.v6.2.symbols.gmt")


## Calculate enrichment as df:s
# DOWN in STAT3
universe_df <- genenames
stat3mt_all_healthy_all_enrich_dn      <- enrich_hypergeometric(df_temp = stat3mt_all_healthy_all_df, dir_temp = "Down")
stat3mt_all_lgl_wtl_enrich_dn          <- enrich_hypergeometric(df_temp = stat3mt_all_lgl_wtl_df, dir_temp = "Down")

stat3mt_cd8_cd48_healthy_all_enrich_dn <- enrich_hypergeometric(stat3mt_cd8_cd48_healthy_all_df, dir_temp = "Down")
stat3mt_cd8_cd48_lgl_wt_enrich_dn      <- enrich_hypergeometric(stat3mt_cd8_cd48_lgl_wt_df, dir_temp = "Down")

stat3mt_cd8_healthy_cd8_enrich_dn      <- enrich_hypergeometric(stat3mt_cd8_healthy_cd8_df, dir_temp = "Down")
stat3mt_cd8_lgl_wt_enrich_dn           <- enrich_hypergeometric(stat3mt_cd8_lgl_wt_df, dir_temp = "Down")

healthy_cd8_healthy_cd4_enrich_dn      <- enrich_hypergeometric(healthy_cd8_healthy_cd4_df, dir_temp = "Down")



enrich_gsea(subset(stat3mt_all_healthy_all_df, direction_cmn == "Down"))
enrich_gsea(subset(stat3mt_all_lgl_wtl_df, direction_cmn == "Down"))

enrich_gsea(subset(stat3mt_cd8_cd48_healthy_all_df, direction_cmn == "Down"))
enrich_gsea(subset(stat3mt_cd8_cd48_lgl_wt_df, direction_cmn == "Down"))

enrich_gsea(subset(stat3mt_cd8_healthy_cd8_df, direction_cmn == "Down"))
enrich_gsea(subset(stat3mt_cd8_lgl_wt_df, direction_cmn == "Down"))

enrich_gsea(subset(healthy_cd8_healthy_cd4_df, direction_cmn == "Down"))





## Write down
write.table(stat3mt_all_healthy_all_enrich_dn, paste0("results/Enrichment/Down_enrich", a_name, ".txt"), sep = "\t")
write.table(stat3mt_all_lgl_wtl_enrich_dn, paste0("results/Enrichment/Down_enrich", b_name, ".txt"), sep = "\t")

write.table(stat3mt_cd8_cd48_healthy_all_enrich_dn, paste0("results/Enrichment/Down_enrich", c_name, ".txt"), sep = "\t")
write.table(stat3mt_cd8_cd48_lgl_wt_enrich_dn, paste0("results/Enrichment/Down_enrich", d_name, ".txt"), sep = "\t")

write.table(stat3mt_cd8_healthy_cd8_enrich_dn, paste0("results/Enrichment/Down_enrich", e_name, ".txt"), sep = "\t")
write.table(stat3mt_cd8_lgl_wt_enrich_dn, paste0("results/Enrichment/Down_enrich", f_name, ".txt"), sep = "\t")

write.table(healthy_cd8_healthy_cd4_enrich_dn, paste0("results/Enrichment/Down_enrich", g_name, ".txt"), sep = "\t")




# UP in STAT3
stat3mt_all_healthy_all_enrich_up      <- enrich_hypergeometric(stat3mt_all_healthy_all_df, dir_temp = "Up")
stat3mt_all_lgl_wtl_enrich_up          <- enrich_hypergeometric(stat3mt_all_lgl_wtl_df, dir_temp = "Up")

stat3mt_cd8_cd48_healthy_all_enrich_up <- enrich_hypergeometric(stat3mt_cd8_cd48_healthy_all_df, dir_temp = "Up")
stat3mt_cd8_cd48_lgl_wt_enrich_up      <- enrich_hypergeometric(stat3mt_cd8_cd48_lgl_wt_df, dir_temp = "Up")

stat3mt_cd8_healthy_cd8_enrich_up      <- enrich_hypergeometric(stat3mt_cd8_healthy_cd8_df, dir_temp = "Up")
stat3mt_cd8_lgl_wt_enrich_up           <- enrich_hypergeometric(healthy_cd8_healthy_cd4_df, dir_temp = "Up")

healthy_cd8_healthy_cd4_enrich_up      <- enrich_hypergeometric(healthy_cd8_healthy_cd4_df, dir_temp = "Up")



enrich_gsea(subset(stat3mt_all_healthy_all_df, direction_cmn == "Up"))
df <- enrich_gsea(subset(stat3mt_all_lgl_wtl_df, direction_cmn == "Up"))

enrich_gsea(subset(stat3mt_cd8_cd48_healthy_all_df, direction_cmn == "Up"))
enrich_gsea(subset(stat3mt_cd8_cd48_lgl_wt_df, direction_cmn == "Up"))

enrich_gsea(subset(stat3mt_cd8_healthy_cd8_df, direction_cmn == "Up"))
enrich_gsea(subset(stat3mt_cd8_lgl_wt_df, direction_cmn == "Up"))

enrich_gsea(subset(healthy_cd8_healthy_cd4_df, direction_cmn == "Up"))


## Write down
write.table(stat3mt_all_healthy_all_enrich_up, paste0("results/Enrichment/Up_enrich", a_name, ".txt"), sep = "\t")
write.table(stat3mt_all_lgl_wtl_enrich_up, paste0("results/Enrichment/Up_enrich", b_name, ".txt"), sep = "\t")

write.table(stat3mt_cd8_cd48_healthy_all_enrich_up, paste0("results/Enrichment/Up_enrich", c_name, ".txt"), sep = "\t")
write.table(stat3mt_cd8_cd48_lgl_wt_enrich_up, paste0("results/Enrichment/Up_enrich", d_name, ".txt"), sep = "\t")

write.table(stat3mt_cd8_healthy_cd8_enrich_up, paste0("results/Enrichment/Up_enrich", e_name, ".txt"), sep = "\t")
write.table(stat3mt_cd8_lgl_wt_enrich_up, paste0("results/Enrichment/Up_enrich", f_name, ".txt"), sep = "\t")

write.table(healthy_cd8_healthy_cd4_enrich_up, paste0("results/Enrichment/Up_enrich", g_name, ".txt"), sep = "\t")






## Plot enrichments as barplots
a <- plot_enrichment(stat3mt_all_healthy_all_enrich_dn, title = paste("Down in",  a_title))
b <- plot_enrichment(stat3mt_all_lgl_wtl_enrich_dn, title = paste("Down in", b_title))

c <- plot_enrichment(stat3mt_cd8_cd48_healthy_all_enrich_dn, title = paste("Down in", c_title))
d <- plot_enrichment(stat3mt_cd8_cd48_lgl_wt_enrich_dn, title = paste("Down in", d_title))

e <- plot_enrichment(stat3mt_cd8_healthy_cd8_enrich_dn, title = paste("Down in", e_title))
f <- plot_enrichment(stat3mt_cd8_lgl_wt_enrich_dn, title = paste("Down in", f_title))

g <- plot_enrichment(healthy_cd8_healthy_cd4_enrich_dn, title = paste("Down in", g_title))


pdf("results/geneset_enrichment_meta_down.pdf", height = 32, width = 8)
grid.arrange(a,b,c,d,e,f,g,ncol=1)
dev.off()


# Save individual plots
grob_list <- list(a,b,c,d,e,f,g)
for(i in 1:7){
  grob_to_pdf(grob_list[[i]], title = "Down_geneset_enrichment", name = name_list[[i]], folder = "Enrichment/", width = 12, height = 12)
}




a <- plot_enrichment(stat3mt_all_healthy_all_enrich_up, title = paste("Up in",  a_title))
b <- plot_enrichment(stat3mt_all_lgl_wtl_enrich_up, title = paste("Up in", b_title))

c <- plot_enrichment(stat3mt_cd8_cd48_healthy_all_enrich_up, title = paste("Up in", c_title))
d <- plot_enrichment(stat3mt_cd8_cd48_lgl_wt_enrich_up, title = paste("Up in", d_title))

e <- plot_enrichment(stat3mt_cd8_healthy_cd8_enrich_up, title = paste("Up in", e_title))
f <- plot_enrichment(stat3mt_cd8_lgl_wt_enrich_up, title = paste("Up in", f_title))

g <- plot_enrichment(healthy_cd8_healthy_cd4_enrich_up, title = paste("Up in", g_title))


pdf("results/geneset_enrichment_meta_Up.pdf", height = 32, width = 8)
grid.arrange(a,b,c,d,e,f,ncol=1)
dev.off()


# Save individual plots
grob_list <- list(a,b,c,d,e,f,g)
for(i in 1:7){
  grob_to_pdf(grob_list[[i]], title = "Up_geneset_enrichment", name = name_list[[i]], folder = "Enrichment/", width = 12, height = 12)
}

