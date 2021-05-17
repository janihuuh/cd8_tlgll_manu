
message("Preprocess TCRab...")

## Read in the sequencing results; the most important files are all_contig and clonotypes -files
contig_files    <- list.files("data/scRNAseq+TCRseq/", recursive = T, pattern = "all_contig_annotations.csv", full.names = T)
clonotype_files <- list.files("data/scRNAseq+TCRseq/", recursive = T, pattern = "clonotypes.csv", full.names = T)

contigs    <- lapply(contig_files, read.delim, sep = ",")
clonotypes <- lapply(clonotype_files, read.delim, sep = ",")

contigs    <- lapply(contig_files, fread)
clonotypes <- lapply(clonotype_files, fread) #, sep = ",")

## Pre-process the data. The function produces two files, barcoded and clonotyped
for(i in 1:length(clonotypes)){

  name <- substr(clonotype_files[i], 1, nchar(clonotype_files[i]) - 15) %>% extractFileName()

  if(!name %in% gsub("_barcoded.txt", "", extractFileName(barcoded_files))){

    message(name)
    preprocess_10X_TCR(contig_file = contigs[[i]], clonotype_file = clonotypes[[i]], prefix = paste0("data/scRNAseq+TCRseq/preprocessed/", name))

   }

}


## Read in the just processed files
clonotyped_files <- list.files("data/scRNAseq+TCRseq/preprocessed/", full.names = T, pattern = "clonotyped")
barcoded_files   <- list.files("data/scRNAseq+TCRseq/preprocessed/", full.names = T, pattern = "barcoded")

clonotyped <- lapply(clonotyped_files, fread)
barcoded   <- lapply(barcoded_files, function(x){fread(x, fill = T) %>% ncol() %>% message(); message(x); fread(x, fill = T) %>% mutate(barcode_uniq = paste0(gsub("_barcoded.txt", "", extractFileName(x)), "_", barcode))}) %>% rbindlist(fill = F)
# barcoded   <- barcoded %>% select(-c(V6:V11))


## Make new clonotype id:s for each of the patient profiled
barcoded$patient <- extractName(barcoded$barcode_uniq)
barcoded_per_pt  <- split(barcoded, f = barcoded$patient)

tot_barcode <- lapply(barcoded_per_pt, newClonotype_df) %>% rbindlist() %>% mutate(new_clonotypes_id = paste0(patient, "_", new_clonotypes_id)) 


tot_barcode$barcode_uniq <- gsub("1189", "fm1189", tot_barcode$barcode_uniq)
tot_barcode$barcode_uniq <- gsub("1890", "fm1890", tot_barcode$barcode_uniq)
tot_barcode$barcode_uniq <- gsub("2311", "fm2311", tot_barcode$barcode_uniq)
tot_barcode$barcode_uniq <- gsub("2477", "fm2477", tot_barcode$barcode_uniq)
tot_barcode$barcode_uniq <- gsub("2481", "fm2481", tot_barcode$barcode_uniq)
tot_barcode$barcode_uniq <- gsub("1338", "fm1338", tot_barcode$barcode_uniq)
tot_barcode$barcode_uniq <- gsub("1656_17", "fm1656_17", tot_barcode$barcode_uniq)
tot_barcode$barcode_uniq <- gsub("1656_18", "fm1656_18", tot_barcode$barcode_uniq)
tot_barcode$barcode_uniq <- gsub("FM1382_2018", "fm 1382", tot_barcode$barcode_uniq)
tot_barcode$barcode_uniq <- gsub("FM805_2015", "fm 805_15", tot_barcode$barcode_uniq)
tot_barcode$barcode_uniq <- gsub("FM805_2018", "fm 805_18", tot_barcode$barcode_uniq)
tot_barcode$barcode_uniq <- substr(tot_barcode$barcode_uniq, 1, nchar(tot_barcode$barcode_uniq) - 2)

## Write down results
write.table(tot_barcode, "data/scRNAseq+TCRseq/preprocessed/tot_barcode.txt",  sep = "\t", row.names = F, quote = F)
tot_barcode[!duplicated(tot_barcode$new_clonotypes_id), ] %>% write.table("data/scRNAseq+TCRseq/preprocessed/tot_uniq_tcrab.txt", sep = "\t", row.names = F, quote = F)






## Dipa edit
barcoded_files <- list.files("data/scRNAseq+TCRseq/preprocessed_dipa/", full.names = T, pattern = "final")
barcoded_final <- lapply(barcoded_files, function(x){fread(x) %>% mutate(name = extractFileName(x) %>% extractName())}) %>% rbindlist() #message(x); fread(x, fill = T) %>% mutate(barcode_uniq = paste0(gsub("_barcoded.txt", "", extractFileName(x)), "_", barcode))}) %>% rbindlist(fill = F)
barcoded_final$barcode_uniq <- paste0(barcoded_final$name, "_", substr(barcoded_final$barcode, 1, nchar(barcoded_final$barcode) - 2))
barcoded_final$Time <- "1"

barcoded_files <- list.files("data/scRNAseq+TCRseq/preprocessed_dipa/", full.names = T, pattern = "join")
barcoded_join <- lapply(barcoded_files, function(x){fread(x) %>% mutate(name = extractFileName(x) %>% extractName())}) %>% rbindlist() #message(x); fread(x, fill = T) %>% mutate(barcode_uniq = paste0(gsub("_barcoded.txt", "", extractFileName(x)), "_", barcode))}) %>% rbindlist(fill = F)
barcoded_join$barcode_uniq <- paste0(barcoded_join$name, "_", barcoded_join$Time, "_", substr(barcoded_join$barcode, 1, nchar(barcoded_join$barcode) - 2))

tot_barcode <- rbind(barcoded_final, barcoded_join)

tot_barcode$barcode_uniq <- gsub("1189", "fm1189", tot_barcode$barcode_uniq)
tot_barcode$barcode_uniq <- gsub("1890", "fm1890", tot_barcode$barcode_uniq)
tot_barcode$barcode_uniq <- gsub("2311", "fm2311", tot_barcode$barcode_uniq)
tot_barcode$barcode_uniq <- gsub("2477", "fm2477", tot_barcode$barcode_uniq)
tot_barcode$barcode_uniq <- gsub("2481", "fm2481", tot_barcode$barcode_uniq)
tot_barcode$barcode_uniq <- gsub("1338", "fm1338", tot_barcode$barcode_uniq)
tot_barcode$barcode_uniq <- gsub("1656_2017", "fm1656_17", tot_barcode$barcode_uniq)
tot_barcode$barcode_uniq <- gsub("1656_2018", "fm1656_18", tot_barcode$barcode_uniq)
tot_barcode$barcode_uniq <- gsub("1382", "fm 1382", tot_barcode$barcode_uniq)
tot_barcode$barcode_uniq <- gsub("805_2015", "fm 805_15", tot_barcode$barcode_uniq)
tot_barcode$barcode_uniq <- gsub("805_2018", "fm 805_18", tot_barcode$barcode_uniq)
# tot_barcode$barcode_uniq <- substr(tot_barcode$barcode_uniq, 1, nchar(tot_barcode$barcode_uniq) - 2)

tot_barcode$Clonotype <- gsub(" ", "", tot_barcode$Clonotype)
tot_barcode$Clonotype <- gsub("C", "c", tot_barcode$Clonotype)
tot_barcode$Clonotype <- paste0(tot_barcode$name, "_", tot_barcode$Clonotype)
tot_barcode$new_clonotypes_id <- tot_barcode$Clonotype

colnames(tot_barcode) <- gsub("\\.x", "_alpha", colnames(tot_barcode))
colnames(tot_barcode) <- gsub("\\.y", "_beta", colnames(tot_barcode))

tot_barcode <- tot_barcode %>% select(clonotype_description, v_gene_alpha:j_gene_alpha, cdr3_alpha, cdr3_nt_alpha, v_gene_beta:j_gene_beta, cdr3_beta, cdr3_nt_beta, name:new_clonotypes_id)
write.table(tot_barcode, "data/scRNAseq+TCRseq/preprocessed/tot_barcode_dipa.txt",  sep = "\t", row.names = F, quote = F)


##### QC
names1 <- substr(colnames(lgl_seurat), 1, nchar(colnames(lgl_seurat)) - 17) %>% unique()
names2 <- substr(tot_barcode$barcode_uniq, 1, nchar(tot_barcode$barcode_uniq) - 19) %>% unique()
names3 <- substr(tot_barcode$barcode_uniq, 1, nchar(tot_barcode$barcode_uniq) - 17) %>% unique()

tot_barcode %>% filter(new_clonotypes_id %in% "FM805_clonotype1") %>% View
lgl_seurat@meta.data %>% filter(new_clonotypes_id %in% "FM805_clonotype1")

barcodes1 <- substr(lgl_seurat$barcode.y, 1, nchar(lgl_seurat$barcode.y) - 2)
barcodes2 <- substr(lgl_seurat$barcode, nchar(lgl_seurat$barcode) - 15, nchar(lgl_seurat$barcode))

table(barcodes2 == barcodes1)

