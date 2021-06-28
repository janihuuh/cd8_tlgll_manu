
getClonalityGenes <- function(object){
  
  clonality_genes <- c(grep("^TRAV", rownames(object), value = T), grep("^TRBV", rownames(object), value = T),
                       grep("^TRGV", rownames(object), value = T), grep("^TRDV", rownames(object), value = T),
                       
                       grep("^TRAJ", rownames(object), value = T), grep("^TRBJ", rownames(object), value = T),
                       grep("^TRGJ", rownames(object), value = T), grep("^TRDJ", rownames(object), value = T),
                       
                       grep("^IGLV", rownames(object), value = T), grep("^IGLC", rownames(object), value = T),
                       grep("^IGLL", rownames(object), value = T), grep("^IGKV", rownames(object), value = T),
                       grep("^IGHV", rownames(object), value = T), grep("^IGKC", rownames(object), value = T),
                       grep("^IGH", rownames(object), value = T),  grep("^IGK", rownames(object), value = T))
  
}



getUnwantedGenes <- function(object){
  
  unwanted_variation <- c(grep("^LINC", rownames(object), value = T), grep("^AC", rownames(object), value = T),
                          grep("^AL", rownames(object), value = T),
                          grep("^MT-", rownames(object), value = T), grep("^RP", rownames(object), value = T))
  
}



## Plot most notable markers
big_markers           <- c("CD3E", "TRAC", "TRBC2",        ## T cell
                           "GNLY", "CD8A", "CD8B", "GZMB", ## CD8+ cell
                           "NKG7", "FCGR3A", "PRF1","KLRG1", "KLRB1", "KLRD1", "NCAM1", ## NK-cell
                           "CD4", "IL7R", "IFNG", "FOXP3", ## CD4+ cell
                           "LYZ", "CD14", "CST3",          ## monocytes
                           "CLEC4C", "PTPRS", "IRF8", "TCF4", ## pDC
                           "MS4A1", "CD19",                   ## b cells
                           "TNFRSF17", "XBP1", "JCHAIN",
                           "MKI67" )                          ## plasma cells
inhibitory_long       <- c("CD28", "ICOS", "CTLA4", "PDCD1", "BTLA", "TNFRSF14", "CD27", "TNFRSF9", "TNFRSF4", "TNFRSF25", "TNFRSF9", "CD2", "CD244", "HAVCR1", "HAVCR2", "TIGIT", "CD226", "CD160", "LAG3", "LAIR1", "CD80")
tcell_activation      <- c("CD69", "CCR7","CD27", "BTLA","CD40LG","IL2RA","CD3E","CD47","EOMES","GNLY", "GZMA","GZMB", "PRF1", "IFNG","CD8A", "CD8B",  "LAMP1","LAG3","CTLA4","HLA-DRA","TNFRSF4","ICOS","TNFRSF9","TNFRSF18")


## Guo et al., Nat Med 2018
naive_markers         <- c("TCF7", "SELL", "LEF1", "CCR7")
cytotoxic_markers     <- c("IL2", "GZMK", "IFNG", "GZMA", "GZMB", "GZMM", "GZMH", "PRF1", "GNLY", "NCAM1")
inhibitory_markers    <- c("LAG3", "TIGIT", "PDCD1", "HAVCR2", "CTLA4")
costimulatory_markers <- c("CD28", "TNFRSF14", "ICOS", "TNFRSF9")
tf_markers            <- c("EOMES", "HOPX",  "TBX21", "ZEB2", "ZNF683", "HIF1A", "ID2", "TOX", "ZBTB16")
treg_markers          <- c("IL2RA", "FOXP3", "LAG3", "IKZF2")
memory_markers        <- c("CD27", "KLRG1", "B3GAT1", "IL2RA")
guo_markers           <- list(naive_markers = naive_markers, cytotoxic_markers = cytotoxic_markers, inhibitory_markers = inhibitory_markers,
                              costimulatory_markers = costimulatory_markers, tf_markers = tf_markers, treg_markers = treg_markers,
                              memory_markers = memory_markers)


## Zhang et al, Nature 2018
cd8_tn                <- c("CCR7", "LEF1", "SELL", "TCF7", "CD27", "CD28", "S1PR1")
cd8_cm                <- c("CCR7", "SELL", "IL7R", "CD27", "CD28", "PRF1", "GZMA", "CCL5", "GPR183", "S1PR1")
cd8_emra              <- c("KLRG1", "CX3CR1", "FCGR3A", "FGFBP2", "PRF1", "GZMH", "TBX21", "EOMES", "S1PR1", "S1PR5")
cd8_em                <- c("GZMK", "CXCR4", "CXCR3", "CD44")
cd8_rm                <- c("CD6", "XCL1", "XCL2", "MYADM", "CAPG", "RORA", "NR4A1", "NR4A2", "NR4A3", "CD69", "ITGAE")
cd8_IEL               <- c("CD160", "KIR2DL4", "TMIGD2", "KLRC1", "KLRC2", "KLRC3", "NR4A1", "NR4A2", "NR4A3", "IKZF2", "ENTPD1", "CD69", "ITGAE")
cd8_tex               <- c("HAVCR2", "PDCD1", "LAYN", "TOX", "IFNG", "GZMB", "MIR155HG", "TNFRSF9", "ITGAE")
cd8_mait              <- c("SLC4A10", "KLRB1", "ZBTB16", "NCR3", "RORC", "RORA")

cd4_tn                <- c("CCR7", "LEF1", "SELL", "TCF7", "CD27", "CD28", "S1PR1")
cd4_tcm_n             <- c("CCR7", "SELL", "PTGER2", "ICAM2", "ANXA1", "ANXA2", "S1PR1")
cd4_emra              <- c("KLRG1", "CX3CR1", "NKG7", "PRF1", "GNLY", "GZMH", "TBX21", "CTSW", "S1PR1", "S1PR5")
cd4_tcm_n             <- c("CCR7", "TCF7", "RGS1", "CD69")
cd4_trm               <- c("CD69", "KLRB1", "PTGER4", "IL7R", "CXCR6",  "NR4A1", "NR4A2", "NR4A3", "MYADM")
cd4_tfh               <- c("CXCR5", "BCL6", "ICA1", "TOX", "TOX2", "IL6ST", "MAGEH1", "BTLA", "ICOS", "PDCD1", "CD200")
cd4_th1               <- c("GZMK", "GZMA", "CCL5", "IFNG", "RUNX3", "EOMES", "CXCR3", "CXCR4", "CD44")
cd4_th17              <- c("IL23R", "RORC", "FURIN", "CTSH", "CCR6", "KLRB1", "CAPG", "ITGAE")
cd4_th1_like          <- c("IFNG", "CXCR3", "BHLHE40", "GZMB", "PDCD1", "HAVCR2", "ICOS", "IGFLR1", "ITGAE")

cd4_treg_blood        <- c("FOXP3", "IL2RA", "IL10RA", "IKZF2", "RTKN2", "CDC25B", "S1PR4")
cd4_tfr               <- c("FOXP3", "IL2RA", "CXCR5", "PDCD1", "IL10", "CCR4", "CD69")
cd4_treg_tum          <- c("FOXP3", "CCR8", "TNFRSF18", "LAYN", "TNFRSF9", "IKZF2", "RTKN2", "CTLA4", "BATF", "IL21R")

zhang_cd8_markers        <- list(cd8_tn, cd8_cm, cd8_emra, cd8_em, cd8_rm, cd8_IEL, cd8_tex, cd8_mait)
names(zhang_cd8_markers) <- c("cd8_tn", "cd8_cm", "cd8_emra", "cd8_em", "cd8_rm", "cd8_IEL", "cd8_tex", "cd8_mait")
zhang_cd8_genes          <- do.call(zhang_cd8_markers, what = "c")

zhang_cd4_markers        <- list(cd4_tn, cd4_tcm_n, cd4_emra, cd4_tcm_n, cd4_trm, cd4_tfh, cd4_th1, cd4_th17, cd4_th1_like, cd4_treg_blood, cd4_treg_tum)
names(zhang_cd4_markers) <- c("cd4_tn", "cd4_tcm_n", "cd4_emra", "cd4_tcm_n", "cd4_trm", "cd4_tfh", "cd4_th1", "cd4_th17", "cd4_th1_like", "cd4_treg_blood", "cd4_treg_tum")
zhang_cd4_genes          <- do.call(zhang_cd4_markers, what = "c")
zhang_markers            <- c(zhang_cd8_markers, zhang_cd4_markers)



#### Van Galen et al Cell 2019

## Undifferntiated
stem_prog             <- c("CD34", "CD38", "KIT", "CD33", "IL3RA", "CD99", "CD96", "PROM1")
promono_mono          <- c("CD68", "CD14", "ITGAM", "ITGAX", "CD64")
hsc                   <- c("MEIS1", "EGR1")
hsc_prog              <- c("MSI2", "CD38", "CD34", "PROM1", "EGFL7")
undiff_mark           <- list(stem_prog = stem_prog, promono_mono = promono_mono, hsc = hsc, hsc_prog = hsc_prog)

## Myeloid
myeloid               <- c("ANPEP", "MPO", "HLA-DRA")
gmp                   <- c("MPO", "ELANE", "CTSG", "AZU1")
promono_mono          <- c("LYST", "LYZ", "CEBPD", "MNDA", "FCER1G")
mono                  <- c("FCN1", "CD14", "C5AR1", "CLEC4A")
cDC                   <- c("CLECL10A", "FCER1A")
pDC                   <- c("CLEC4C", "PTPRS", "IRF8", "TCF4")
myeloid_mark          <- list(myeloid = myeloid, gmp = gmp, promono_mono = promono_mono, mono = mono, cDC = cDC, pDC = pDC)

## Erythroid
early_late_ery        <- c("CSF1", "KIT", "HBB", "HBD")
late_ery              <- "GYPA"
erythroid_mark        <- list(early_late_ery = early_late_ery, late_ery = late_ery)

## Lymphoid
pro_b                 <- c("CD24", "EBF1", "MME", "VPREB1")
pro_b_b               <- c("PAX5", "CD19", "CD79")
b                     <- c("MS4A1", "BANK1")
plasma                <- c("MZB1", "IGLL5", "JCHAIN")

t_ctl                 <- c("CD3D", "CD3G")
t_ctl_nk              <- c("IL32", "IL7R")
t_cell                <- c("TCF7")
ctl_nk                <- c("CCL5")
ctl                   <- c("GZMK")

nk                    <- c("KLRB1", "KLRD1", "GZMB", "NCAM1")
lymphoid_mark         <- list(pro_b = pro_b, pro_b_b = pro_b_b, b = b, plasma = plasma, t_ctl = t_ctl, t_ctl_nk = t_ctl_nk, t_cell = t_cell, ctl_nk = ctl_nk, ctl = ctl, nk = nk)

van_galen_markers     <- c(undiff_mark, myeloid_mark, erythroid_mark, lymphoid_mark)


# Pfefferle et al. biorxiv 2019
is_nk                 <- c("CD3E", "KLRB1", "KLRD1", "GZMB", "NCAM1")
nk_cd56_bright        <- c("NCAM1", "KLRC1", "MYC", "BACH2")
nk_cd56_dim           <- c("FCERG3A", "PRDM1", "ZEB2", "MAF")
nk_adaptive           <- c("KLRC2", "KLF2")
nk_mark               <- list(is_nk = is_nk, nk_cd56_bright = nk_cd56_bright, nk_cd56_dim = nk_cd56_dim, nk_adaptive = nk_adaptive)





# ======

marker_list           <- list(big_markers, naive_markers, cytotoxic_markers, inhibitory_markers, costimulatory_markers, tf_markers, treg_markers, memory_markers) # , th1_markers, th2_markers, th17_markers, th1_flow_markers, th2_flow_markers, th17_flow_markers, thF_flow_markers)
marker_genes          <- do.call("c", marker_list)
marker_amount         <- do.call("c", lapply(marker_list, length))

marker_list           <- list(big_markers, naive_markers, cytotoxic_markers, costimulatory_markers, tf_markers, treg_markers, inhibitory_long, tcell_activation)
marker_genes          <- do.call("c", marker_list)
marker_amount         <- do.call("c", lapply(marker_list, length))
split_known           <- rep(c("Overall", "Naive", "Cytotoxic", "Costimulatory", "TF", "Treg", "Inhibitory", "T_cell_activation"), marker_amount)
