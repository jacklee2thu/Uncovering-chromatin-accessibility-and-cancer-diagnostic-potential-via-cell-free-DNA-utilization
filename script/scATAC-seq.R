##ArchR
options(stringsAsFactors=F)
# setwd('E:/科研/singlecell_ATAC_case')
setwd('/Share2/home/lanxun3/jacklee/cfDNA/singlecell_ATAC/temp_data')
inputFiles<-as.character(c("/Share2/home/lanxun3/jacklee/cfDNA/singlecell_ATAC/Stomach_normal1/GSM5589404_stomach_SM-CHLWL_rep1_fragments_final.tsv.gz",
"/Share2/home/lanxun3/jacklee/cfDNA/singlecell_ATAC/Stomach_normal2/GSM5589405_stomach_SM-IOBHV_rep1_fragments_final.tsv.gz",
"/Share2/home/lanxun3/jacklee/cfDNA/singlecell_ATAC/Stomach_normal3/GSM5589406_stomach_SM-JF1NP_rep1_fragments_final.tsv.gz",
"/Share2/home/lanxun3/jacklee/cfDNA/singlecell_ATAC/Stomach_normal4/GSM5589407_stomach_SM-JF1O3_rep1_fragments_final.tsv.gz",
"/Share2/home/lanxun3/jacklee/cfDNA/singlecell_ATAC/GIST/GSM5276954_3E4D1L_ATAC_fragments.tsv.gz"))
names(inputFiles)<-c("GSM5589404","GSM5589405","GSM5589406","GSM5589407","GIST")
library(Cairo)
library(ArchR)
rhdf5::h5disableFileLocking()
set.seed(1)
addArchRThreads(threads = 24) 
addArchRGenome("hg38")
ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  sampleNames = names(inputFiles),
  minTSS = 4, #Dont set this too high because you can always increase later
  minFrags = 1000, 
  addTileMat = TRUE,
  addGeneScoreMat = TRUE
)
ArrowFiles<-as.character(c("GSM5589404.arrow","GSM5589405.arrow","GSM5589406.arrow","GSM5589407.arrow","GIST.arrow"))
doubScores <- addDoubletScores(
  input = ArrowFiles,
  k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
  knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search.
  LSIMethod = 1
)
proj <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = "singlecellATAC",
  copyArrows = F #This is recommened so that you maintain an unaltered copy for later usage.
)
proj <- filterDoublets(ArchRProj = proj)##去除双细胞
###降维与聚类/PeakMatrix/TileMatrix
proj <- addIterativeLSI(ArchRProj = proj, useMatrix = "TileMatrix", 
name = "IterativeLSI",clusterParams = list(resolution = c(20), sampleCells = 10000, maxClusters = 25, n.start= 10))
proj_addHarmony <- addHarmony(
    ArchRProj = proj,
    reducedDims = "IterativeLSI",
    name = "Harmony",
    groupBy = "Sample"
)
proj_addClusters <- addClusters(input = proj_addHarmony, reducedDims = "IterativeLSI",maxClusters =50,knnAssign = 25)
proj_addClusters <- addClusters(input = proj, reducedDims = "IterativeLSI",maxClusters =50)
##画图
proj_addUMAP <- addUMAP(ArchRProj = proj_addClusters, reducedDims = "IterativeLSI")
proj_umap <- saveArchRProject(ArchRProj = proj_addUMAP)
proj_addUMAP <- loadArchRProject(path = "/Share2/home/lanxun3/jacklee/cfDNA/singlecell_ATAC/temp_data/singlecellATAC")

p1 <- plotEmbedding(ArchRProj = proj_addUMAP, colorBy = "cellColData", name = "Sample", embedding = "UMAP",
pal=c('#00BEC2','#3E4F9F','#E21A19','#006400','#DAA520','#B22222'))
p2 <- plotEmbedding(ArchRProj = proj_addUMAP, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
ggAlignPlots(p1, type = "h")
plotPDF(p1, name = "singlecellATAC-UMAP-sample.pdf",
        ArchRProj = proj_addUMAP, addDOC = FALSE, width = 5, height = 5)
####sample type
temp_type<-proj_addUMAP@cellColData@rownames
temp_type<-unlist(lapply(strsplit(temp_type,"#"),function(x){return(x[[1]])}))
temp_type[temp_type=='GIST']<-'Tumor'
temp_type[temp_type=='GSM5589404']<-'Normal'
temp_type[temp_type=='GSM5589405']<-'Normal'
temp_type[temp_type=='GSM5589406']<-'Normal'
temp_type[temp_type=='GSM5589407']<-'Normal'
proj_addUMAP@cellColData@listData$sampletype<-temp_type
p4 <- plotEmbedding(ArchRProj = proj_addUMAP, colorBy = "cellColData", name = "sampletype", embedding = "UMAP",
pal=c('#447CAA','#841B1F'))
ggAlignPlots(p4, type = "h")
plotPDF(p4, name = "singlecellATAC-UMAP-sampletype.pdf",
        ArchRProj = proj_addUMAP, addDOC = FALSE, width = 5, height = 5)
		
#cell type
proj_addUMAP@cellColData@listData$celltype<-proj_addUMAP@cellColData@listData$Clusters
proj_addUMAP@cellColData@listData$celltype[proj_addUMAP@cellColData@listData$celltype=='C1']<-'DC'
proj_addUMAP@cellColData@listData$celltype[proj_addUMAP@cellColData@listData$celltype=='C2']<-'normal Epithelial'
proj_addUMAP@cellColData@listData$celltype[proj_addUMAP@cellColData@listData$celltype=='C3']<-'normal Epithelial'
proj_addUMAP@cellColData@listData$celltype[proj_addUMAP@cellColData@listData$celltype=='C4']<-'normal Epithelial'
proj_addUMAP@cellColData@listData$celltype[proj_addUMAP@cellColData@listData$celltype=='C5']<-'mast'
proj_addUMAP@cellColData@listData$celltype[proj_addUMAP@cellColData@listData$celltype=='C6']<-'T'
proj_addUMAP@cellColData@listData$celltype[proj_addUMAP@cellColData@listData$celltype=='C7']<-'B'
proj_addUMAP@cellColData@listData$celltype[proj_addUMAP@cellColData@listData$celltype=='C8']<-'Endocrine'
proj_addUMAP@cellColData@listData$celltype[proj_addUMAP@cellColData@listData$celltype=='C9']<-'fibroblast'
proj_addUMAP@cellColData@listData$celltype[proj_addUMAP@cellColData@listData$celltype=='C10']<-'Mac'
proj_addUMAP@cellColData@listData$celltype[proj_addUMAP@cellColData@listData$celltype=='C11']<-'B'
proj_addUMAP@cellColData@listData$celltype[proj_addUMAP@cellColData@listData$celltype=='C12']<-'T'
proj_addUMAP@cellColData@listData$celltype[proj_addUMAP@cellColData@listData$celltype=='C13']<-'T'
proj_addUMAP@cellColData@listData$celltype[proj_addUMAP@cellColData@listData$celltype=='C14']<-'T'
proj_addUMAP@cellColData@listData$celltype[proj_addUMAP@cellColData@listData$celltype=='C15']<-'T'
proj_addUMAP@cellColData@listData$celltype[proj_addUMAP@cellColData@listData$celltype=='C16']<-'Mac'
proj_addUMAP@cellColData@listData$celltype[proj_addUMAP@cellColData@listData$celltype=='C17']<-'fibroblast'
proj_addUMAP@cellColData@listData$celltype[proj_addUMAP@cellColData@listData$celltype=='C18']<-'Tumor Epithelial'
p3 <- plotEmbedding(ArchRProj = proj_addUMAP, colorBy = "cellColData", name = "celltype", embedding = "UMAP",
pal=c('B'='#1AB4B8','T'='#ED877F','Endocrine'='#E06AA4','normal Epithelial'='#EEDD82','mast'='#8585B9','Tumor Epithelial'='#DDB056',
'DC'='#9E9F20','Mac'="#7DC05B",'fibroblast'="#7DC6AC"))
ggAlignPlots(p3, type = "h")
plotPDF(p3, name = "singlecellATAC-UMAP-cell.pdf",
        ArchRProj = proj_addUMAP, addDOC = FALSE, width = 5, height = 5)
#######识别marker gene
markerGenes  <- c('KRT18','EPCAM','CD3D','CD68','AIF1','CHGA','CDH5',
'ENG','PECAM1','DCN','ID2','MS4A1','CD79A','IGLL5','CPA3')

library(presto)
markersGS <- getMarkerFeatures(
    ArchRProj = proj_addUMAP, 
    useMatrix = "GeneScoreMatrix", 
    groupBy = "Clusters",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon"
)

index<-match(markerGenes,markersGS@elementMetadata@listData$name)

FC_marker<-markersGS@assays@data@listData$Log2FC[index,]
p_marker<-markersGS@assays@data@listData$Pval[index,]
rownames(FC_marker)<-markerGenes
colnames(FC_marker)<-paste('C',c(1:18))
rownames(p_marker)<-markerGenes
colnames(p_marker)<-paste('C',c(1:19))

cell_meta<-data.frame(proj_addUMAP@cellColData@listData$celltype,proj_addUMAP@cellColData@rownames)
colnames(cell_meta)<-c('cell_type','cell_barcode')
save(cell_meta,file='cell_meta.Rdata')

###
markersGS <- getMarkerFeatures(
    ArchRProj = proj, 
    useMatrix = "GeneScoreMatrix", 
    groupBy = "celltype",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon"
)

index<-match(markerGenes,markersGS@elementMetadata@listData$name)

FC_marker<-markersGS@assays@data@listData$Log2FC[index,]
rownames(FC_marker)<-markerGenes
colnames(FC_marker)<-names(table(proj@cellColData@listData$celltype))

FC_marker[FC_marker<= -2.744572]<- -2.744572
library(pheatmap)
pdf('marker_genes.pdf')
pheatmap(FC_marker,cluster_cols = F,cluster_rows = F,color = colorRampPalette(c("#236DAD","white","#690D19"))(50),
show_colnames = T,show_rownames = T,main="marker_gene")
dev.off()

markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.01 & Log2FC >= 1")


###########展示marker gene 热图
heatmapGS <- plotMarkerHeatmap(
  seMarker = markersGS, 
  cutOff = "FDR <= 0.01 & Log2FC >= 1", 
  labelMarkers = markerGenes,
  transpose = TRUE
)
ComplexHeatmap::draw(heatmapGS, heatmap_legend_side = "bot", annotation_legend_side = "bot")

plotPDF(heatmapGS, name = "GeneScores-Marker-Heatmap", width = 8, height = 6, ArchRProj =proj, addDOC = FALSE)


##给cluster增加gene 得分
proj_addImputeWeights <- addImputeWeights(proj_addUMAP)
markerGenes  <- c('KRT18','EPCAM','MUC6','CD3D','CD68','AIF1','CHGA','CDH5','ENG','PECAM1','DCN','ID2','MS4A1','CD79A','IGLL5','CPA3')
markerGenes  <- c('RAE1','GPR85','CALU','AK9','TRAF2','TENT2','ABCA10','LRRIQ1','SPECC1','MEIS2')
final_signature_index<-c('BLOC1S3','TRAF2','CALU','ABCA10','ZNF775','NDUFV2','OAZ2','RAE1','C17orf53','LRRIQ1','SPECC1','TMEM203','WDR34','AK9',
'MFSD2A','TCOF1','CD47','C2CD5','FRS2','SYTL2','MYO9A','TCF4','NPPC','TMTC2','TMX3','MEIS2','PABPC3','GMFB','GCA','ATXN2L','AATF',
'C3orf67','RNF138','OLA1','PPA2','BDKRB1','ARF4','SENP7','PRR15','GPR85')
proj_umap <- saveArchRProject(ArchRProj = proj_addUMAP)
p <- plotEmbedding(
    ArchRProj = proj_addImputeWeights, 
    colorBy = "GeneScoreMatrix", 
    name = markerGenes, 
    embedding = "UMAP",
	pal=c('#DBF0F5','#BEBEBE','#FF0000'),
    imputeWeights = getImputeWeights(proj_addImputeWeights), log2Norm = FALSE
)
#Rearrange for grid plotting
p2 <- lapply(p, function(x){
    x + guides(color = FALSE, fill = FALSE) + 
    theme_ArchR(baseSize = 6.5) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    theme(
        axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(), 
        axis.text.y=element_blank(), 
        axis.ticks.y=element_blank()
    )
})
do.call(cowplot::plot_grid, c(list(ncol = 3),p2))
plotPDF(plotList = p, 
    name = "cfDNA-Marker-Genes-W-Imputation.pdf", 
    ArchRProj = proj_addImputeWeights, 
    addDOC = FALSE, width = 5, height = 5)

ArchRProj = proj_addImputeWeights
colorBy = "GeneScoreMatrix"
name = markerGenes 
imputeWeights = getImputeWeights(proj_addImputeWeights)

colorMat <- ArchR:::.getMatrixValues(ArchRProj = ArchRProj, name = name,matrixName = colorBy, log2Norm = FALSE)
colorMat <- imputeMatrix(mat = as.matrix(colorMat),imputeWeights = imputeWeights)
save(colorMat,file='colorMat.Rdata')
	
#########gene score 在Epithelial和Epithelial_n中比较
setwd('/Share2/home/lanxun3/jacklee/cfDNA/singlecell_ATAC/temp_data')
load('colorMat.Rdata')
markerGenes  <- c('RAE1','GPR85','CALU','AK9','TRAF2','TENT2','ABCA10','LRRIQ1','SPECC1','MEIS2')

load('cell_meta.Rdata')

marker_TSS<-colorMat
Duc_cell<-cell_meta[grep('Tumor',cell_meta$cell_type),]$cell_barcode
Duc_cell<-intersect(Duc_cell,colnames(colorMat))
Duc_n_cell<-cell_meta[grep('normal',cell_meta$cell_type),]$cell_barcode
Duc_n_cell<-intersect(Duc_n_cell,colnames(colorMat))
marker_TSS_p<-apply(marker_TSS,1,function(x){a<-wilcox.test(x[Duc_cell],x[Duc_n_cell]);return(a$p.value)})
names(marker_TSS_p)<-rownames(marker_TSS)
#####Violin plot
Gene_score<-c()
for(i in 1:dim(marker_TSS)[1]){
temp_score<-c(marker_TSS[i,Duc_cell],marker_TSS[i,Duc_n_cell])
Gene_score<-c(Gene_score,temp_score)
}
df_fc<-data.frame(Genescore=Gene_score,gene=factor(c(rep(names(marker_TSS_p),each=length(temp_score))),levels=names(marker_TSS_p)),
type=factor(c(rep(c(rep("tumor",length(Duc_cell)),rep("normal",length(Duc_n_cell))),10)),levels=c('tumor','normal')))
df_fc<-df_fc[df_fc$gene%in%names(marker_TSS_p)[6:10],]
df_fc$combine_name<-factor(paste(df_fc$gene,df_fc$type),levels=c("RAE1 tumor","RAE1 normal","GPR85 tumor","GPR85 normal","CALU tumor",
"CALU normal","AK9 tumor","AK9 normal","TRAF2 tumor","TRAF2 normal",
"TENT2 tumor","TENT2 normal","ABCA10 tumor","ABCA10 normal","LRRIQ1 tumor","LRRIQ1 normal","SPECC1 tumor","SPECC1 normal","MEIS2 tumor","MEIS2 normal"))



library(ggplot2)
library(ggpubr)
color = c("#BE4E4D","#35A1D3")
pdf('STAD_gene_score2.pdf',width=20,height=10)
my_comparisons <- list(c("tumor", "normal"))
df_fc %>% ggplot(aes(x=combine_name,y=Genescore,fill=type))+
  geom_boxplot(width=0.5,position=position_dodge(width=1),outlier.colour = NA)+
  theme_bw()+
  theme(legend.position="none",
        axis.text.x = element_text(hjust = 0.5, vjust = 0.5),
        axis.text.y = element_text(hjust = 0.5, vjust = 0.5),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text( size=rel(1)),
        panel.border = element_blank(),
        axis.line = element_line(colour = "black",size=1)
  )+scale_fill_manual(values = color)+
  stat_compare_means( comparisons = my_comparisons)+
  stat_compare_means(label.y = 0.8)
dev.off()




##########从browser track里看基因开放性
p <- plotBrowserTrack(
    ArchRProj = proj_addImputeWeights, 
    groupBy = "celltype", 
    geneSymbol = markerGenes, 
    upstream = 30000,
    downstream = 30000
)
grid::grid.newpage()
grid::grid.draw(p$CD14)
plotPDF(plotList = p, 
    name = "Plot-Tracks-Marker-Genes.pdf", 
    ArchRProj = proj_addImputeWeights, 
    addDOC = FALSE, width = 5, height = 5)
