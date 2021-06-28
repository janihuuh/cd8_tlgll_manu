
library(SCENIC)
library(AUCell)

lgl_clones_seurat <- readRDS("results/lgl_clones_seurat2.rds")

unwanted_genes <- getUnwantedGenes(lgl_clones_seurat)
clonality_genes <- getClonalityGenes(lgl_clones_seurat)

exprMat <- lgl_clones_seurat@assays$RNA@counts
df <- exprMat %>% as.matrix()

### Initialize settings
org <- "hgnc" # or hgnc, or dmel
dbDir <- "cisTarget_databases" # RcisTarget databases location
myDatasetTitle <- "SCENIC on T-LGLL" # choose a name for your analysis
data(defaultDbNames)
dbs <- defaultDbNames[[org]]
scenicOptions <- initializeScenic(org=org, dbDir=dbDir, dbs=dbs, datasetTitle=myDatasetTitle, nCores=10)
scenicOptions <- readRDS(file="int/scenicOptions.Rds") 


### Co-expression network
genesKept <- geneFiltering(df, scenicOptions)
genesKept <- genesKept[!genesKept %in% c(unwanted_genes, clonality_genes)]
exprMat_filtered <- df[genesKept, ]
exprMat_log <- log2(df+1)
exprMat_filtered_log <- log2(exprMat_filtered+1) 

runGenie3(exprMat_filtered_log, scenicOptions, resumePreviousRun = T)
# runCorrelation(exprMat_filtered, scenicOptions)
runCorrelation(exprMat_filtered_log, scenicOptions)

### Build and score the GRN
scenicOptions@settings$dbs <- scenicOptions@settings$dbs
scenicOptions <- runSCENIC_1_coexNetwork2modules(scenicOptions)
scenicOptions <- runSCENIC_2_createRegulons(scenicOptions) #, coexMethod=c("top5perTarget")) # Toy run settings
scenicOptions <- runSCENIC_3_scoreCells(scenicOptions, exprMat_log)

# Optional: Binarize activity
# aucellApp <- plotTsne_AUCellApp(scenicOptions, exprMat_log)
# savedSelections <- shiny::runApp(aucellApp)
# newThresholds <- savedSelections$thresholds
# scenicOptions@fileNames$int["aucell_thresholds",1] <- "int/newThresholds.Rds"
# saveRDS(newThresholds, file=getIntName(scenicOptions, "aucell_thresholds"))
scenicOptions <- runSCENIC_4_aucell_binarize(scenicOptions)
tsneAUC(scenicOptions, aucType="AUC") # choose settings

# Export:
saveRDS(cellInfo, file=getDatasetInfo(scenicOptions, "cellInfo")) # Temporary, to add to loom
export2loom(scenicOptions, exprMat)

# To save the current status, or any changes in settings, save the object again:
saveRDS(scenicOptions, file="int/scenicOptions.Rds") 

### Exploring output 
# Check files in folder 'output'
# Browse the output .loom file @ http://scope.aertslab.org

# output/Step2_MotifEnrichment_preview.html in detail/subset:
motifEnrichment_selfMotifs_wGenes <- loadInt(scenicOptions, "motifEnrichment_selfMotifs_wGenes")
tableSubset <- motifEnrichment_selfMotifs_wGenes[highlightedTFs=="Sox8"]
viewMotifs(tableSubset) 

# output/Step2_regulonTargetsInfo.tsv in detail: 
regulonTargetsInfo <- loadInt(scenicOptions, "regulonTargetsInfo")
tableSubset <- regulonTargetsInfo[TF=="Stat6" & highConfAnnot==TRUE]
viewMotifs(tableSubset) 

regulonTargetsInfo %>% filter(spearCor > 0.5)

# Cell-type specific regulators (RSS): 
regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
rss <- calcRSS(AUC=getAUC(regulonAUC), cellAnnotation=cellInfo[colnames(regulonAUC), "CellType"], )
rssPlot <- plotRSS(rss)
plotly::ggplotly(rssPlot$plot)


scenicOptions@fileNames




## Visualize results

lgl_clone_seurat
umap.coord <- lgl_clone_seurat@reductions$latent_umap@cell.embeddings
colnames(umap.coord) <- c("UMAP_1", "UMAP_2")

tfs <- c("STAT3",
         "JUNB",
         "JUN",
         
         "FOS",
         "NFKB1",
         "NFKB2")

dir.create("results/scenic/")

pdf("results/scenic/umap_tfs.pdf", width = 9, height = 12)
par(mfrow = c(4,3))
AUCell_plotTSNE(tSNE=umap.coord, exprMat=exprMat, cellsAUC=selectRegulons(regulonAUC, tfs), plots = "binaryAUC", cex = 0.25, exprCols = c("red3")) #, thresholds = T, offColor = "gray95")
dev.off()


tfs <- c("FOS (187g)", "JUN (106g)", "IRF1 (155g)", "JUNB (15g)", "STAT3 (12g)")
tfs <- c("FOS", "JUN", "IRF1", "JUNB", "STAT3")

pdf("results/scenic/umap_tfs2.pdf", width = 5, height = 12)
par(mfrow = c(5,2))
AUCell_plotTSNE(tSNE=umap.coord, exprMat=exprMat, cellsAUC=selectRegulons(regulonAUC, tfs), plots = c("binaryAUC", "AUC"), cex = 0.25, exprCols = c("red3")) #, thresholds = T, offColor = "gray95")
dev.off()


png("results/scenic/umap_tfs2.png", width = 5, height = 12, units = "in", res = 1e3)
par(mfrow = c(5,2))
AUCell_plotTSNE(tSNE=umap.coord, exprMat=exprMat, cellsAUC=selectRegulons(regulonAUC, tfs), plots = c("binaryAUC", "AUC"), cex = 0.25, exprCols = c("red3")) #, thresholds = T, offColor = "gray95")
dev.off()




cellInfo <- data.frame(CellType=Idents(lgl_clones_seurat))

regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]

regulonActivity_byCellType <- sapply(split(rownames(cellInfo), cellInfo$CellType),
                                     function(cells) rowMeans(getAUC(regulonAUC)[,cells]))

regulonActivity_byCellType_Scaled <- t(scale(t(regulonActivity_byCellType), center = T, scale=T))

# ComplexHeatmap::Heatmap(regulonActivity_byCellType, name="Regulon activity")
ComplexHeatmap::Heatmap(regulonActivity_byCellType_Scaled, name="Regulon activity")

regulonActivity_byCellType_Scaled2 <- t(scale(t(regulonActivity_byCellType[,1:5]), center = T, scale=T))

regulonActivity_byCellType_Scaled2 <- t(scale(t(regulonActivity_byCellType[rownames(regulonActivity_byCellType) %in% topRegulators$Regulon,1:5]), center = T, scale=T))

regulonActivity_byCellType_Scaled2 <- t(scale(t(regulonActivity_byCellType[rownames(regulonActivity_byCellType) %in% topRegulators$Regulon,1:5]), center = T, scale=T))
regulonActivity_byCellType_Scaled2 <- t(scale(t(regulonActivity_byCellType[rownames(regulonActivity_byCellType) %in% grep("extended", rownames(regulonActivity_byCellType), value = T, invert = T),1:5]), center = T, scale=T))

set.seed(123)
pdf("results/scenic/heatmap_km4.pdf", width = 6, height = 15)
ComplexHeatmap::Heatmap(regulonActivity_byCellType_Scaled2, name="Regulon activity", col = c("dodgerblue", "white", "red3"), border = T, km = 4)
dev.off()



topRegulators <- reshape2::melt(regulonActivity_byCellType_Scaled2)
colnames(topRegulators) <- c("Regulon", "CellType", "RelativeActivity")
topRegulators <- topRegulators[which(topRegulators$RelativeActivity>0),]
length(topRegulators$Regulon)








minPerc <- .7
binaryRegulonActivity <- loadInt(scenicOptions, "aucell_binary_nonDupl")
cellInfo_binarizedCells <- cellInfo[which(rownames(cellInfo)%in% colnames(binaryRegulonActivity)),, drop=FALSE]
regulonActivity_byCellType_Binarized <- sapply(split(rownames(cellInfo_binarizedCells), cellInfo_binarizedCells$CellType), 
                                               function(cells) rowMeans(binaryRegulonActivity[,cells, drop=FALSE]))
binaryActPerc_subset <- regulonActivity_byCellType_Binarized[which(rowSums(regulonActivity_byCellType_Binarized>minPerc)>0),]

ComplexHeatmap::Heatmap(binaryActPerc_subset[,1:5], name="Regulon activity (%)", col = c("white","pink","red"))

set.seed(123)
pdf("results/scenic/heatmap_km4_bin.pdf", width = 6, height = 12)
ComplexHeatmap::Heatmap(binaryActPerc_subset[,1:5], name="Regulon\nactivity (%)", col = c("dodgerblue", "white", "red3"), border = T, km = 4)
dev.off()

set.seed(123)
pdf("results/scenic/heatmap_km4_bin2.pdf", width = 4.5, height = 6)
ComplexHeatmap::Heatmap(binaryActPerc_subset[rownames(binaryActPerc_subset) %in% grep("extended", rownames(binaryActPerc_subset), value = T, invert = T), 1:5], name="Regulon\nactivity (%)", col = c("dodgerblue", "white", "red3"), border = T, km = 4)
dev.off()


     