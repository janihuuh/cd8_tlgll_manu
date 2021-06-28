
#### TCR analysis of LGL patient Cohort ######

## pre processing of TCR data doen for each sample individually
tcr_file <- read.csv2("data/tcrseq/all_contig_annotations2015.csv")
# basic filtering
tcr_file <- subset(tcr_file, is_cell==TRUE & high_confidence==TRUE & full_length==TRUE & raw_clonotype_id != "None" & cdr3 != "None" & cdr3_nt != "None" & chain != "Multi")
## there are 9237 cells with TCR information
cell_ids <- unique(tcr_file$barcode)
length(cell_ids)
### there are 6905 cells with unique barcodes

## save all cells that appear more than twice with TCR data
count_of_barcodes <- count(tcr_file,barcode)
max(count_of_barcodes$n)
multi_barcode_cells <- subset(count_of_barcodes, n>2)
multi_barcodes <- multi_barcode_cells$barcode
MULTI_BARCODE_TCR_data <- subset(tcr_file, barcode %in% multi_barcodes)
write.csv2(MULTI_BARCODE_TCR_data, "lgl_analysis/TCR_analysis/805/2015/MULTI_barcode_cells.csv")

######remove all cells with multi barcodes ###
tcr_file <- subset(tcr_file, !(barcode %in% multi_barcodes))
cell_ids <- unique(tcr_file$barcode)
length(cell_ids)
###Find all cells that have both alpha and beta chains
alpha_tcr <- subset(tcr_file,chain == "TRA" )
beta_tcr <- subset(tcr_file,chain == "TRB" )
intersect_AB_barcodes <- intersect(subset(tcr_file,chain == "TRA" )$barcode,subset(tcr_file,chain == "TRB" )$barcode)
## 3470 cells have paired TCR information #####

### Merging
merged_meta_clonotype <- merge(alpha_tcr, beta_tcr, by.x = "barcode", by.y = "barcode")
identical(merged_meta_clonotype$raw_clonotype_id.x, merged_meta_clonotype$raw_clonotype_id.y)
### Merged data####
write.csv2(merged_meta_clonotype, "lgl_analysis/TCR_analysis/805/2015/AB_merged.csv",row.names = FALSE )

tcr_ab_merged <- read.csv2("lgl_analysis/TCR_analysis/805/2015/AB_merged.csv")
max(count(tcr_ab_merged,cdr3.x)$n)
max(count(tcr_ab_merged,cdr3.y)$n)
max(count(tcr_ab_merged,cdr3_nt.x)$n)
max(count(tcr_ab_merged,cdr3_nt.y)$n)




### adding new conotype annotation depending on amino acid seq ####
###############################################################
tcr <-unite(tcr_ab_merged, "clonotype_description", c("cdr3.x", "cdr3.y"), sep = "_", remove = FALSE)
max(count(tcr,clonotype_description)$n)
count_clonotypes <- count(tcr,vars=clonotype_description)
count_clonotypes <- count_clonotypes [order(count_clonotypes$n, decreasing = TRUE),]
nrow(count_clonotypes)
length(count_clonotypes)
count_clonotypes [,3] <- 0
colnames(count_clonotypes) <- c("clonotype_description","Count","Clonotype")
for (i in 1:nrow(count_clonotypes)) {
  count_clonotypes [i,3] <- paste("Clonotype",i)
}
final_merge <- merge(tcr, count_clonotypes, all.x = TRUE, all.y = TRUE)
write_csv2(final_merge,"lgl_analysis/TCR_analysis/805/2015/final_tcr_file.csv")
write_csv2(count_clonotypes,"lgl_analysis/TCR_analysis/805/2015/final_count_file.csv")

##################################################################
##################################################################


##################################################################
##################################################################
################### Annotating cells with Leukemic and non leukemic clones  ###########################
tcr1338 <- read.csv2("lgl_analysis/TCR_analysis/1338/final_tcr_file.csv")
clones <- c("Clonotype 1")
tcr1338 <- subset(tcr1338, Clonotype %in% clones)
unique(tcr1338$Clonotype)
tcr1338 <- as.data.table(tcr1338)
tcr1338 <- tcr1338[,barcode := gsub("\\-1","",barcode)]
tcr1338 <- tcr1338[,barcode :=paste0("fm1338_", barcode)]
tcr1338$Clone_Description <- tcr1338[,38]
tcr1338 <- tcr1338[,Clone_Description := gsub("Clonotype 1","sure_lgll",Clone_Description)]
tcr1338$sample <- "FM1338"
rm(clones)
##############################################

##############################################
tcr1382 <- read.csv2("lgl_analysis/TCR_analysis/1382/final_tcr_file.csv")
clones <- c("Clonotype 1", "Clonotype 2")
tcr1382 <- subset(tcr1382, Clonotype %in% clones)
unique(tcr1382$Clonotype)
tcr1382 <- as.data.table(tcr1382)
tcr1382 <- tcr1382[,barcode := gsub("\\-1","",barcode)]
tcr1382 <- tcr1382[,barcode :=paste0("fm 1382_", barcode)]
tcr1382$Clone_Description <- tcr1382[,38]
tcr1382 <- tcr1382[,Clone_Description := gsub("Clonotype 1","unsure_lgll",Clone_Description)]
tcr1382 <- tcr1382[,Clone_Description := gsub("Clonotype 2","sure_lgll",Clone_Description)]
tcr1382$sample <- "FM1382"
rm(clones)
##############################################

##############################################
tcr1189 <- read.csv2("lgl_analysis/TCR_analysis/1189/final_tcr_file.csv")
clones <- c("Clonotype 1", "Clonotype 2","Clonotype 3")
tcr1189 <- subset(tcr1189, Clonotype %in% clones)
unique(tcr1189$Clonotype)
tcr1189 <- as.data.table(tcr1189)
tcr1189 <- tcr1189[,barcode := gsub("\\-1","",barcode)]
tcr1189 <- tcr1189[,barcode :=paste0("fm1189_", barcode)]
tcr1189$Clone_Description <- tcr1189[,38]
tcr1189 <- tcr1189[,Clone_Description := gsub("Clonotype 1","unsure_lgll",Clone_Description)]
tcr1189 <- tcr1189[,Clone_Description := gsub("Clonotype 2","unsure_lgll",Clone_Description)]
tcr1189 <- tcr1189[,Clone_Description := gsub("Clonotype 3","unsure_lgll",Clone_Description)]
tcr1189$sample <- "FM1189"
rm(clones)
##############################################

##############################################
tcr1890 <- read.csv2("lgl_analysis/TCR_analysis/1890/final_tcr_file.csv")
clones <- c("Clonotype 1", "Clonotype 2")
tcr1890 <- subset(tcr1890, Clonotype %in% clones)
unique(tcr1890$Clonotype)
tcr1890 <- as.data.table(tcr1890)
tcr1890 <- tcr1890[,barcode := gsub("\\-1","",barcode)]
tcr1890 <- tcr1890[,barcode :=paste0("fm1890_", barcode)]
tcr1890$Clone_Description <- tcr1890[,38]
tcr1890 <- tcr1890[,Clone_Description := gsub("Clonotype 1","sure_lgll",Clone_Description)]
tcr1890 <- tcr1890[,Clone_Description := gsub("Clonotype 2","unsure_lgll",Clone_Description)]
tcr1890$sample <- "FM1890"
rm(clones)
##############################################

##############################################
tcr805_15 <- read.csv2("lgl_analysis/TCR_analysis/805/2015/final_tcr_file.csv")
clones <- c("Clonotype 1", "Clonotype 2")
tcr805_15 <- subset(tcr805_15, Clonotype %in% clones)
unique(tcr805_15$Clonotype)
tcr805_15 <- as.data.table(tcr805_15)
tcr805_15 <- tcr805_15[,barcode := gsub("\\-1","",barcode)]
tcr805_15 <- tcr805_15[,barcode :=paste0("fm 805_15_", barcode)]
tcr805_15$Clone_Description <- tcr805_15[,38]
tcr805_15 <- tcr805_15[,Clone_Description := gsub("Clonotype 1","sure_lgll",Clone_Description)]
tcr805_15 <- tcr805_15[,Clone_Description := gsub("Clonotype 2","unsure_lgll",Clone_Description)]
tcr805_15$sample <- "FM805_15"
rm(clones)
##############################################

##############################################
tcr805_18 <- read.csv2("lgl_analysis/TCR_analysis/805/2018/final_tcr_file.csv")
clones <- c("Clonotype 1", "Clonotype 2","Clonotype 4")
tcr805_18 <- subset(tcr805_18, Clonotype %in% clones)
unique(tcr805_18$Clonotype)
tcr805_18 <- as.data.table(tcr805_18)
tcr805_18 <- tcr805_18[,barcode := gsub("\\-1","",barcode)]
tcr805_18 <- tcr805_18[,barcode :=paste0("fm 805_18_", barcode)]
tcr805_18$Clone_Description <- tcr805_18[,38]
tcr805_18 <- tcr805_18[,Clone_Description := gsub("Clonotype 1","unsure_lgll",Clone_Description)]
tcr805_18 <- tcr805_18[,Clone_Description := gsub("Clonotype 2","unsure_lgll",Clone_Description)]
tcr805_18 <- tcr805_18[,Clone_Description := gsub("Clonotype 4","sure_lgll",Clone_Description)]
tcr805_18$sample <- "FM805_18"
rm(clones)
##############################################

##############################################
tcr1656_17 <- read.csv2("lgl_analysis/TCR_analysis/1656/2017/final_tcr_file.csv")
clones <- c("Clonotype 1")
tcr1656_17 <- subset(tcr1656_17, Clonotype %in% clones)
unique(tcr1656_17$Clonotype)
tcr1656_17 <- as.data.table(tcr1656_17)
tcr1656_17 <- tcr1656_17[,barcode := gsub("\\-1","",barcode)]
tcr1656_17 <- tcr1656_17[,barcode :=paste0("fm1656_17_", barcode)]
tcr1656_17$Clone_Description <- tcr1656_17[,38]
tcr1656_17 <- tcr1656_17[,Clone_Description := gsub("Clonotype 1","sure_lgll",Clone_Description)]
tcr1656_17$sample <- "FM1656_17"
rm(clones)
##############################################

##############################################
tcr1656_18 <- read.csv2("lgl_analysis/TCR_analysis/1656/2018/final_tcr_file.csv")
clones <- c("Clonotype 1")
tcr1656_18 <- subset(tcr1656_18, Clonotype %in% clones)
unique(tcr1656_18$Clonotype)
tcr1656_18 <- as.data.table(tcr1656_18)
tcr1656_18 <- tcr1656_18[,barcode := gsub("\\-1","",barcode)]
tcr1656_18 <- tcr1656_18[,barcode :=paste0("fm1656_18_", barcode)]
tcr1656_18$Clone_Description <- tcr1656_18[,38]
tcr1656_18 <- tcr1656_18[,Clone_Description := gsub("Clonotype 1","sure_lgll",Clone_Description)]
tcr1656_18$sample <- "FM1656_18"
rm(clones)
##############################################

##############################################
tcr2481 <- read.csv2("lgl_analysis/TCR_analysis/2481/final_tcr_file.csv")
clones <- c("Clonotype 1","Clonotype 2")
tcr2481 <- subset(tcr2481, Clonotype %in% clones)
unique(tcr2481$Clonotype)
tcr2481 <- as.data.table(tcr2481)
tcr2481 <- tcr2481[,barcode := gsub("\\-1","",barcode)]
tcr2481 <- tcr2481[,barcode :=paste0("fm2481_", barcode)]
tcr2481$Clone_Description <- tcr2481[,38]
tcr2481 <- tcr2481[,Clone_Description := gsub("Clonotype 1","sure_lgll",Clone_Description)]
tcr2481 <- tcr2481[,Clone_Description := gsub("Clonotype 2","unsure_lgll",Clone_Description)]
tcr2481$sample <- "FM2481"
rm(clones)
##############################################

##############################################
tcr2477 <- read.csv2("lgl_analysis/TCR_analysis/2477/final_tcr_file.csv")
clones <- c("Clonotype 1","Clonotype 2")
tcr2477 <- subset(tcr2477, Clonotype %in% clones)
unique(tcr2477$Clonotype)
tcr2477 <- as.data.table(tcr2477)
tcr2477 <- tcr2477[,barcode := gsub("\\-1","",barcode)]
tcr2477 <- tcr2477[,barcode :=paste0("fm2477_", barcode)]
tcr2477$Clone_Description <- tcr2477[,38]
tcr2477 <- tcr2477[,Clone_Description := gsub("Clonotype 1","sure_lgll",Clone_Description)]
tcr2477 <- tcr2477[,Clone_Description := gsub("Clonotype 2","unsure_lgll",Clone_Description)]
tcr2477$sample <- "FM2477"
rm(clones)
##############################################

##############################################
tcr2311 <- read.csv2("lgl_analysis/TCR_analysis/2311/final_tcr_file.csv")
clones <- c("Clonotype 1","Clonotype 2")
tcr2311 <- subset(tcr2311, Clonotype %in% clones)
unique(tcr2311$Clonotype)
tcr2311 <- as.data.table(tcr2311)
tcr2311 <- tcr2311[,barcode := gsub("\\-1","",barcode)]
tcr2311 <- tcr2311[,barcode :=paste0("fm2311_", barcode)]
tcr2311$Clone_Description <- tcr2311[,38]
tcr2311 <- tcr2311[,Clone_Description := gsub("Clonotype 1","sure_lgll",Clone_Description)]
tcr2311 <- tcr2311[,Clone_Description := gsub("Clonotype 2","sure_lgll",Clone_Description)]
tcr2311$sample <- "FM2311"
rm(clones)

##############################################

all_clones <- rbind(tcr1189,tcr1338,tcr1382,tcr1656_17,tcr1656_18,tcr1890,tcr2311,tcr2477,tcr2481,tcr805_15,tcr805_18)
write.csv2(all_clones, "lgl_analysis/re_annotated_clones/all_clones.csv", row.names = FALSE)
sure_lglls <- subset(all_clones, Clone_Description == "sure_lgll")
unsure_lglls <- subset(all_clones, Clone_Description == "unsure_lgll")
write.csv2(sure_lglls, "lgl_analysis/re_annotated_clones/sure_lgll.csv", row.names = FALSE)
write.csv2(unsure_lglls, "lgl_analysis/re_annotated_clones/unsure_lgll.csv", row.names = FALSE)


################## Adding augmented B chain clones ###
augment2481 <- read.csv2("lgl_analysis/TCR_analysis/2481/augmentbeta_cnltp1_2.csv.csv")
augment2481 <- as.data.table(augment2481)
augment2481 <- augment2481[,barcode := gsub("\\-1","",barcode)]
augment2481 <- augment2481[,barcode :=paste0("fm2481_", barcode)]
augment2481$sample <- "FM2481"

augment1656_17 <- read.csv2("lgl_analysis/TCR_analysis/1656/2017/Beta_TCR_Clonotype1.csv")
augment1656_17 <- as.data.table(augment1656_17)
augment1656_17 <- augment1656_17[,barcode := gsub("\\-1","",barcode)]
augment1656_17 <- augment1656_17[,barcode :=paste0("fm1656_17_", barcode)]
augment1656_17$sample <- "FM1656_17"

augment1656_18 <- read.csv2("lgl_analysis/TCR_analysis/1656/2018/Beta_TCR_Clonotype1.csv")
augment1656_18 <- as.data.table(augment1656_18)
augment1656_18 <- augment1656_18[,barcode := gsub("\\-1","",barcode)]
augment1656_18 <- augment1656_18[,barcode :=paste0("fm1656_18_", barcode)]
augment1656_18$sample <- "FM1656_18"

augmented_allclones <- rbind(augment1656_17,augment1656_18,augment2481)
write.csv2(augmented_allclones, "lgl_analysis/re_annotated_clones/augment_clones.csv", row.names = FALSE)

beta_sure_lglls <- subset(augmented_allclones, Clone_description == "sure_lgll")
beta_unsure_lglls <- subset(augmented_allclones, Clone_description == "unsure_lgll")
write.csv2(beta_sure_lglls, "lgl_analysis/re_annotated_clones/beta_sure_lgll.csv", row.names = FALSE)
write.csv2(beta_unsure_lglls, "lgl_analysis/re_annotated_clones/beta_unsure_lgll.csv", row.names = FALSE)
