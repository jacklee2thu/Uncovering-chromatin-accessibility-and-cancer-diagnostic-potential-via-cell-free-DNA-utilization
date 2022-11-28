####比对
cd /Share2/home/lanxun3/jacklee/cfDNA/PDAC_ATAC/ATAC_1836591

/Share2/home/zhaolab/xiap/software/cellranger-atac-1.2.0/cellranger-atac count --id=ATAC_1836591 \
			--reference=/Share2/home/zhaolab/xiap/refdata/refdata-cellranger-atac-hg19-1.2.0  \
			--fastqs=/Share2/home/zhaolab/data/PDAC_ATAC/atac-seq_2_1836591_FKDN202559715-1A \
			--sample=FKDL202607050-1a-SI-NA-C2 \
			--localcores=50 \
			--localmem=50
			
cd /Share2/home/lanxun3/jacklee/cfDNA/singlecell_ATAC/Stomach_normal4
gunzip -c GSM5589407_stomach_SM-JF1O3_rep1_fragments.bed.gz > GSM5589407_stomach_SM-JF1O3_rep1_fragments.bed
########normal pancreas scATAC
wget https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5589391&format=file&file=GSM5589391%5Fpancreas%5FSM%2DADRUQ%5Frep1%5Ffragments%2Ebed%2Egz
cd /home/lijie/SPOT_ATAC/GSM5589391_pancreas_SM-ADRUQ_rep1
gunzip GSM5589391_pancreas_SM-ADRUQ_rep1_fragments.bed.gz
####liftOver hg38tohg19
cd /home/lijie/SPOT_ATAC/GSM5589391_pancreas_SM-ADRUQ_rep1
/home/lijie/soft_ware/liftOver/liftOver GSM5589394_pancreas_SM-JF1O6_rep1_fragments.bed /home/lijie/soft_ware/liftOver/hg38ToHg19.over.chain.gz \
GSM5589394_pancreas_SM-JF1O6_rep1_fragments_hg19.bed hg19_unmap.bed

########
cd /home/lijie/SPOT_ATAC/GSM5589394_pancreas_SM-JF1O6_rep1
sort -k 1,1 -k2,2n GSM5589394_pancreas_SM-JF1O6_rep1_fragments_hg19.tsv > GSM5589394_pancreas_SM-JF1O6_rep1_fragments_final.tsv
bgzip -c GSM5589394_pancreas_SM-JF1O6_rep1_fragments_final.tsv > GSM5589394_pancreas_SM-JF1O6_rep1_fragments_final.tsv.gz
tabix -b 2 -e 3 -p bed GSM5589394_pancreas_SM-JF1O6_rep1_fragments_final.tsv.gz
rm *.tsv

#######TSSEnrichment, CreateRegionPileupMatrix,计算细胞的TSS score
#######采用TSS upstream 1k和downstream 1k[normalize c(1:100, 1902:2001)]来计算cfDNA TSS2K对应的TSS score(-log2transform)
#######采用TSS upstream 150和downstream 50[normalize c(1:100, 1902:2001)]来计算cfDNA TSS2K对应的TSS score(-log2transform)
##
##第二种方法
#######采用reads(TSS upstream 1k和downstream 1k)/2k / [normalize reads(c(1:2000, 4000:6000))/4k]来计算cfDNA TSS2K对应的TSS score
#######采用reads(TSS upstream 150和downstream 50)/200 / [normalize reads(c(1:2000, 4000:6000))/4k]来计算cfDNA TSS2K对应的TSS score

#########createArrowFiles, 
#########识别细胞类型,ArchR
options(stringsAsFactors=F)
setwd('/home/lijie/SPOT_ATAC/fragments')
inputFiles<-as.character(c("/home/lijie/SPOT_ATAC/ATAC_1836591/fragments.tsv.gz",
"/home/lijie/SPOT_ATAC/ATAC_T001835227/fragments.tsv.gz",
"/home/lijie/SPOT_ATAC/ATAC_T001837519/fragments.tsv.gz",
"/home/lijie/SPOT_ATAC/ATAC_T001837519_n/fragments.tsv.gz",
"/home/lijie/SPOT_ATAC/GSM5589391_pancreas_SM-ADRUQ_rep1/GSM5589391_pancreas_SM-ADRUQ_rep1_fragments_final.tsv.gz",
"/home/lijie/SPOT_ATAC/GSM5589392_pancreas_SM-IOBHS_rep1/GSM5589392_pancreas_SM-IOBHS_rep1_fragments_final.tsv.gz",
"/home/lijie/SPOT_ATAC/GSM5589393_pancreas_SM-JF1NS_rep1/GSM5589393_pancreas_SM-JF1NS_rep1_fragments_final.tsv.gz",
"/home/lijie/SPOT_ATAC/GSM5589394_pancreas_SM-JF1O6_rep1/GSM5589394_pancreas_SM-JF1O6_rep1_fragments_final.tsv.gz"))
names(inputFiles)<-c("ATAC1","ATAC2",'ATAC3','ATAC_n','ATAC_p1','ATAC_p2','ATAC_p3','ATAC_p4')

setwd('/Share2/home/lanxun3/jacklee/cfDNA/PDAC_ATAC/fragments')
inputFiles<-as.character(c("/Share2/home/lanxun3/jacklee/cfDNA/PDAC_ATAC/ATAC_1836591/ATAC_1836591/outs/fragments.tsv.gz",
"/Share2/home/lanxun3/jacklee/cfDNA/PDAC_ATAC/ATAC_T001835227/ATAC_T001835227/outs/fragments.tsv.gz",
"/Share2/home/lanxun3/jacklee/cfDNA/PDAC_ATAC/ATAC_T001837519/ATAC_T001837519_t/outs/fragments.tsv.gz",
"/Share2/home/lanxun3/jacklee/cfDNA/PDAC_ATAC/ATAC_T001837519/ATAC_T001837519_n/outs/fragments.tsv.gz"))
names(inputFiles)<-c("ATAC1","ATAC2",'ATAC3','ATAC_n')
library(Cairo)
library(ArchR)
rhdf5::h5disableFileLocking()
set.seed(1)
addArchRThreads(threads = 23) 
addArchRGenome("hg19")
ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  sampleNames = names(inputFiles),
  minTSS = 4, #Dont set this too high because you can always increase later
  minFrags = 1000, 
  addTileMat = TRUE,
  addGeneScoreMat = TRUE
)
ArrowFiles<-as.character(c('ATAC1.arrow','ATAC2.arrow','ATAC3.arrow','ATAC_n.arrow'))
# ArrowFiles<-as.character(c('ATAC1.arrow','ATAC2.arrow','ATAC3.arrow','ATAC_n.arrow','ATAC_p1.arrow','ATAC_p2.arrow','ATAC_p3.arrow','ATAC_p4.arrow'))
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

proj_umap <- saveArchRProject(ArchRProj = proj)
proj <- loadArchRProject(path = "/home/lijie/SPOT_ATAC/fragments/singlecellATAC")

###降维与聚类
proj <- addIterativeLSI(ArchRProj = proj, useMatrix = "TileMatrix",name = "IterativeLSI",varFeatures = 30000,clusterParams = list(resolution = c(10), sampleCells = 10000, maxClusters = 20, n.start= 10),force = TRUE)
proj <- addHarmony(
    ArchRProj = proj,
    reducedDims = "IterativeLSI",
    name = "Harmony",
    groupBy = "Sample",
	force = TRUE
)
proj <- addClusters(input = proj, reducedDims = "IterativeLSI",maxClusters = 50,force = TRUE)

##画图
proj <- addUMAP(ArchRProj = proj, reducedDims = "IterativeLSI",force = TRUE)
proj_umap <- saveArchRProject(ArchRProj = proj)
proj <- loadArchRProject(path = "/home/lijie/SPOT_ATAC/fragments/singlecellATAC")

p1 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Sample", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
ggAlignPlots(p2, type = "h")
plotPDF(p1, name = "Plot-UMAP-Sample.pdf",
        ArchRProj = proj, addDOC = FALSE, width = 5, height = 5)
		
#########
################sample type
temp_type<-proj@cellColData@rownames
temp_type<-unlist(lapply(strsplit(temp_type,"#"),function(x){return(x[[1]])}))
temp_type[temp_type=='ATAC1']<-'Tumor'
temp_type[temp_type=='ATAC2']<-'Tumor'
temp_type[temp_type=='ATAC3']<-'Tumor'
temp_type[temp_type=='ATAC_n']<-'Normal'
temp_type[temp_type=='ATAC_p1']<-'Normal'
temp_type[temp_type=='ATAC_p2']<-'Normal'
temp_type[temp_type=='ATAC_p3']<-'Normal'
temp_type[temp_type=='ATAC_p4']<-'Normal'
proj@cellColData@listData$sampletype<-temp_type
p4 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "sampletype", embedding = "UMAP",
pal=c('Normal'='#447CAA','Tumor'='#841B1F'))
ggAlignPlots(p4, type = "h")
plotPDF(p4, name = "PDAC-UMAP-sampletype.pdf",
        ArchRProj = proj, addDOC = FALSE, width = 5, height = 5)


#cell type
proj@cellColData@listData$celltype<-proj@cellColData@listData$Clusters
proj@cellColData@listData$celltype[proj@cellColData@listData$celltype=='C1']<-'B'
proj@cellColData@listData$celltype[proj@cellColData@listData$celltype=='C2']<-'B'
proj@cellColData@listData$celltype[proj@cellColData@listData$celltype=='C3']<-'B'
proj@cellColData@listData$celltype[proj@cellColData@listData$celltype=='C4']<-'T'
proj@cellColData@listData$celltype[proj@cellColData@listData$celltype=='C5']<-'T'
proj@cellColData@listData$celltype[proj@cellColData@listData$celltype=='C6']<-'NKT'
proj@cellColData@listData$celltype[proj@cellColData@listData$celltype=='C7']<-'T'
proj@cellColData@listData$celltype[proj@cellColData@listData$celltype=='C8']<-'ENDOCRINE'
proj@cellColData@listData$celltype[proj@cellColData@listData$celltype=='C9']<-'DUC_n'
proj@cellColData@listData$celltype[proj@cellColData@listData$celltype=='C10']<-'ACINAR'
proj@cellColData@listData$celltype[proj@cellColData@listData$celltype=='C11']<-'ACINAR'
proj@cellColData@listData$celltype[proj@cellColData@listData$celltype=='C12']<-'DUC'
proj@cellColData@listData$celltype[proj@cellColData@listData$celltype=='C13']<-'DUC'
proj@cellColData@listData$celltype[proj@cellColData@listData$celltype=='C14']<-'DUC_n'
proj@cellColData@listData$celltype[proj@cellColData@listData$celltype=='C15']<-'DC'
proj@cellColData@listData$celltype[proj@cellColData@listData$celltype=='C16']<-'MAC'
proj@cellColData@listData$celltype[proj@cellColData@listData$celltype=='C17']<-'FIBROBLAST'
proj@cellColData@listData$celltype[proj@cellColData@listData$celltype=='C18']<-'ENDOTHELIAL'
proj@cellColData@listData$celltype[proj@cellColData@listData$celltype=='C19']<-'FIBROBLAST'


p3 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "celltype", embedding = "UMAP",
pal=c('B'='#1AB4B8','T'='#ED877F','NKT'='#E06AA4','ENDOCRINE'='#E06AA4','DUC_n'='#EE82EE','ACINAR'='#87CEFA','DUC'='#A020F0',
'DC'='#9E9F20','MAC'="#7DC05B",'FIBROBLAST'="#7DC6AC",'ENDOTHELIAL'="#F4A460"))
ggAlignPlots(p3, type = "h")
plotPDF(p3, name = "Plot-UMAP-cell.pdf",
        ArchRProj = proj, addDOC = FALSE, width = 5, height = 5)

proj_umap <- saveArchRProject(ArchRProj = proj)
proj <- loadArchRProject(path = "/Share2/home/lanxun3/jacklee/cfDNA/PDAC_ATAC/fragments/singlecellATAC")
#######识别marker gene
markerGenes  <- c('CD3D','CD4','CD8A','CD8B','TCF7','CCR7','FOXP3','CTLA4','PDCD1','GZMK','GZMB','CXCL13',
'FCGR3A','TYROBP','CD19','MS4A1','CD79A','CD79B','IGLL5','CD27','MKI67','CDH5','DCN','LUM','CHGB','PRSS1','KRT19','AIF1','ITGAX','CD1C','CD68','MPO')
markerGenes  <- c('CD3D','FCGR3A','CD19','MS4A1','CD79A','CD79B','IGLL5','KRT19','MKI67','CDH5','DCN','LUM','CHGB','PRSS1','EPCAM','AIF1','ITGAX','CD1C','CD68')

library(presto)
markersGS <- getMarkerFeatures(
    ArchRProj = proj, 
    useMatrix = "GeneScoreMatrix", 
    groupBy = "Clusters",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon"
)

index<-match(markerGenes,markersGS@elementMetadata@listData$name)

FC_markers<-markersGS@assays[index,]

FC_marker<-FC_markers@data[[1]]
p_marker<-FC_markers@data[[4]]
mean_marker<-FC_markers@data[[2]]
rownames(FC_marker)<-markerGenes
colnames(FC_marker)<-paste('C',c(1:19))
rownames(p_marker)<-markerGenes
colnames(p_marker)<-paste('C',c(1:19))
rownames(mean_marker)<-markerGenes
colnames(mean_marker)<-paste('C',c(1:19))

cell_meta<-data.frame(proj@cellColData@listData$celltype,proj@cellColData@rownames)
colnames(cell_meta)<-c('cell_type','cell_barcode')
save(cell_meta,file='cell_meta.Rdata')
##ATAC_p2有10000+细胞，缩减到7000个
setwd('/home/lijie/SPOT_ATAC/fragments')
load('cell_meta.Rdata')
options(stringsAsFactors=F)
GSM5589392_pancreas_SM_IOBHS_rep1_cellmatrix<-cell_meta[grep('ATAC_p2',cell_meta$cell_barcode),]
set.seed(1234)
sample_index<-sample(GSM5589392_pancreas_SM_IOBHS_rep1_cellmatrix$cell_barcode,size = 7000)
new_cell_map<-cell_meta[-grep('ATAC_p2',cell_meta$cell_barcode),]
GSM5589392_pancreas_SM_IOBHS_rep1_cellmatrixn<-GSM5589392_pancreas_SM_IOBHS_rep1_cellmatrix[match(sample_index,GSM5589392_pancreas_SM_IOBHS_rep1_cellmatrix$cell_barcode),]
cell_meta<-rbind(new_cell_map,GSM5589392_pancreas_SM_IOBHS_rep1_cellmatrixn)
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

FC_markers<-markersGS@assays[index,]
mean_marker<-FC_markers@data[[2]]
rownames(mean_marker)<-markerGenes
colnames(mean_marker)<-names(table(proj@cellColData@listData$celltype))


library(pheatmap)
pdf('marker_genes.pdf')
pheatmap(t(scale(t(mean_marker))),cluster_cols = T,cluster_rows = F,color = colorRampPalette(c("#236DAD","white","#690D19"))(50),
show_colnames = T,show_rownames = T,main="marker_gene")
dev.off()

markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.01 & Log2FC >= 1.5")
save(markerList,file='markerList_celltype.Rdata')###########各个细胞类型的marker
cell_type_marker<-unique(c(markerList$ACINAR$name,markerList$B$name,markerList$DC$name,markerList$DUC$name,markerList$DUC_n$name,markerList$ENDOCRINE$name,
markerList$ENDOTHELIAL$name,markerList$FIBROBLAST$name,markerList$MAC$name,markerList$NKT$name,markerList$T$name))

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
proj <- addImputeWeights(proj)
markerGenes  <- c('AURKC','PRSS33','ZFPM2','GAP43','CFTR','ZNF616','PRDM6','BMF','ZNF580','BORCS7')
markerGenes  <- c("AURKC","PRSS33","ZFPM2","GAP43","CFTR","ZNF616","PRDM6","BMF","ZNF580","BRMS1","PARP8","NFIA","BORCS7","CELF4","NLRP4",
"MDFI","CCHCR1","ZNF667","MYO10","ZFPL1","MAP3K19","ENSA","B3GNT3","GREB1L","ZNF165","CYP2U1","VRK1","PIFO","MYBPC3","DRG2","COL11A1","GGA1",
"TRIM44","TPSAB1","ADH1A","ZFAND2B","RAB11B","CTF1","CAPN10","ARMC9","CNOT2","GSDMB","SYVN1","ZSCAN21","AKAP6","DOCK5")
markerGenes  <- c('CD3D','FCGR3A','CD19','MS4A1','CD79A','CD79B','IGLL5','KRT19','MKI67','CDH5','DCN','LUM','CHGB','PRSS1','EPCAM','AIF1','ITGAX','CD1C','CD68')
markerGenes  <- c('CFTR','ZNF616','NFIA','MDFI','ZNF667','MYO10','ZFPL1','ENSA','PIFO','GGA1')
p <- plotEmbedding(
    ArchRProj = proj, 
    colorBy = "GeneScoreMatrix", 
    name = markerGenes, 
    embedding = "UMAP",
	pal=c('#DBF0F5','#BEBEBE','#FF0000'),
    imputeWeights = getImputeWeights(proj), log2Norm = FALSE
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
    name = "Plot-UMAP-PDAC_markergenes-W-Imputation.pdf", 
    ArchRProj = proj, 
    addDOC = FALSE, width = 5, height = 5)

ArchRProj = proj
colorBy = "GeneScoreMatrix"
name = markerGenes 
imputeWeights = getImputeWeights(proj)

colorMat <- ArchR:::.getMatrixValues(ArchRProj = ArchRProj, name = name,matrixName = colorBy, log2Norm = FALSE)
colorMat <- imputeMatrix(mat = as.matrix(colorMat),imputeWeights = imputeWeights)
save(colorMat,file='/home/lijie/SPOT_ATAC/SPOT_test/colorMat.Rdata')

##########从browser track里看基因开放性
p <- plotBrowserTrack(
    ArchRProj = proj, 
    groupBy = "celltype", 
    geneSymbol = markerGenes, 
    upstream = 50000,
    downstream = 50000
)
grid::grid.newpage()
grid::grid.draw(p$CD14)
plotPDF(plotList = p, 
    name = "Plot-Tracks-Marker-Genes.pdf", 
    ArchRProj = proj, 
    addDOC = FALSE, width = 5, height = 5)

#########gene score 在DUC和DUC_n中比较
setwd('/home/lijie/SPOT_ATAC/SPOT_test')
load('/home/lijie/SPOT_ATAC/SPOT_test/colorMat.Rdata')
markerGenes  <- c('AURKC','PRSS33','ZFPM2','GAP43','CFTR','ZNF616','PRDM6','BMF','ZNF580','C10orf32')
markerGenes  <- c('CAPN10','ZFAND2B','VRK1','RAB11B','BORCS7','DRG2','ZNF165','ARMC9','DOCK5','SYVN1','CYP2U1','ZSCAN21','ZFPL1','CCHCR1','CNOT2',
'TRIM44','MDFI','TPSAB1','ZNF667','ADH1A','ENSA','BMF','NFIA','GSDMB','NLRP4','BRMS1','GGA1','AKAP6','AURKC','PARP8','PIFO','ZNF580','CELF4',
'PRSS33','CTF1','MYBPC3','COL11A1','B3GNT3','ZFPM2','MYO10','CFTR','MAP3K19','GAP43','GREB1L','ZNF616','PRDM6')
markerGenes  <- c('CFTR','ZNF616','NFIA','MDFI','ZNF667','MYO10','ZFPL1','ENSA','PIFO','GGA1')
load('/home/lijie/SPOT_ATAC/SPOT_test/cell_meta.Rdata')

marker_TSS<-colorMat
Duc_cell<-cell_meta[grep('DUC',cell_meta$cell_type),]$cell_barcode
Duc_cell<-intersect(Duc_cell,colnames(colorMat))
Duc_n_cell<-cell_meta[grep('DUC_n',cell_meta$cell_type),]$cell_barcode
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
df_fc$combine_name<-factor(paste(df_fc$gene,df_fc$type),levels=c("CFTR tumor","CFTR normal","ZNF616 tumor","ZNF616 normal","NFIA tumor",
"NFIA normal","MDFI tumor","MDFI normal","ZNF667 tumor","ZNF667 normal",
"MYO10 tumor","MYO10 normal","ZFPL1 tumor","ZFPL1 normal","ENSA tumor","ENSA normal","PIFO tumor","PIFO normal","GGA1 tumor","GGA1 normal"))



library(ggplot2)
library(ggpubr)
color = c("#BE4E4D","#35A1D3")
pdf('PDAC_gene_score2.pdf',width=20,height=10)
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
  stat_compare_means(label.y = 4)
dev.off()

########排除ACINAR cell,然后计算markercelltype list
options(stringsAsFactors=F)
setwd('/home/lijie/SPOT_ATAC/fragments')
library(Cairo)
library(ArchR)
rhdf5::h5disableFileLocking()
set.seed(1)
addArchRThreads(threads = 23) 
addArchRGenome("hg19")
proj <- loadArchRProject(path = "/home/lijie/SPOT_ATAC/fragments/singlecellATAC")
load('/home/lijie/SPOT_ATAC/SPOT_test/cell_meta.Rdata')
ACINAR_index<-grep('ACINAR',cell_meta$cell_type)###删除ACINAR cell
cell_meta<-cell_meta[-ACINAR_index,]
new_sample<-cell_meta$cell_barcode
pro_new<-subsetArchRProject(
  ArchRProj = proj,
  cells = new_sample,
  outputDirectory = "/home/lijie/SPOT_ATAC/fragments/new",
  dropCells = TRUE,
  logFile = NULL,
  threads = getArchRThreads(),
  force = TRUE
)

# pro_new <- addIterativeLSI(ArchRProj = pro_new, useMatrix = "TileMatrix",force = TRUE, name = "IterativeLSI")
# pro_new <- addHarmony(
    # ArchRProj = pro_new,
    # reducedDims = "IterativeLSI",
    # name = "Harmony",
    # groupBy = "Sample",
	# force = TRUE
# )
# pro_new <- addClusters(input = pro_new,force = TRUE, reducedDims = "IterativeLSI")
# pro_new <- addUMAP(ArchRProj = pro_new,force = TRUE, reducedDims = "IterativeLSI")
setwd('/home/lijie/SPOT_ATAC/fragments/new')
pro_new_umap <- saveArchRProject(ArchRProj = pro_new)
pro_new_umap <- loadArchRProject(path = "/home/lijie/SPOT_ATAC/fragments/new")
p1 <- plotEmbedding(ArchRProj = pro_new_umap, colorBy = "cellColData", name = "Sample", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj =pro_new_umap, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
ggAlignPlots(p1, type = "h")
plotPDF(p1, name = "Plot-UMAP-Sample.pdf",
        ArchRProj = proj, addDOC = FALSE, width = 5, height = 5)
###
markersGS <- getMarkerFeatures(
    ArchRProj = pro_new_umap, 
    useMatrix = "GeneScoreMatrix", 
    groupBy = "celltype",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon"
)


markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.01 & Log2FC >= 2")
# markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.01 & Log2FC >= 2.5")
save(markerList,file='markerList_celltype.Rdata')###########各个细胞类型的marker
# save(markerList,file='markerList_celltype_2.5.Rdata')###########各个细胞类型的marker
cell_type_marker<-unique(c(markerList$B$name,markerList$DC$name,markerList$DUC$name,markerList$DUC_n$name,markerList$ENDOCRINE$name,
markerList$ENDOTHELIAL$name,markerList$FIBROBLAST$name,markerList$MAC$name,markerList$NKT$name,markerList$T$name))
length(cell_type_marker)


##########
#############Signac提取TSSEnrichment
options(stringsAsFactors=F)
library(Seurat)
library(Signac)
library(Matrix)
library(ggplot2)
library(patchwork)
library(GenomicRanges)
set.seed(1234)
#####制作GRanges数据
load('/home/lijie/SPOT_ATAC/SPOT_test/tss_pro_table.Rdata')
tss_pro_table$strand[tss_pro_table$strand==1]<-'+'
tss_pro_table$strand[tss_pro_table$strand== -1]<-'-'
annotations <- GRanges(seqnames = Rle(tss_pro_table$chromosome),
ranges = IRanges(start=tss_pro_table$tss, end = tss_pro_table$tss),
strand = Rle(tss_pro_table$strand),
gene_name = tss_pro_table$gene_name,
gene_biotype = tss_pro_table$gene_type,
gene_id = tss_pro_table$ensemble_id)

tss.positions <- annotations
tss.positions <- Extend(x = tss.positions, upstream = 3000, 
        downstream = 3000, from.midpoint = TRUE)

#########提取细胞的TSS region reads
load('/home/lijie/SPOT_ATAC/SPOT_test/cell_meta.Rdata')
GSM5589392_pancreas_SM_IOBHS_rep1_cellmatrix<-cell_meta[grep('ATAC_p2',cell_meta$cell_barcode),]
GSM5589392_pancreas_SM_IOBHS_rep1_cell<-unlist(lapply(strsplit(GSM5589392_pancreas_SM_IOBHS_rep1_cellmatrix$cell_barcode,'#'),function(x){return(x[[2]])}))

setwd('/home/lijie/SPOT_ATAC/GSM5589392_pancreas_SM-IOBHS_rep1')
cellmap<-GSM5589392_pancreas_SM_IOBHS_rep1_cell
names(cellmap)<-cellmap
region<-tss.positions
tabix.file <- Rsamtools:::TabixFile(file = '/home/lijie/SPOT_ATAC/GSM5589392_pancreas_SM-IOBHS_rep1/GSM5589392_pancreas_SM-IOBHS_rep1_fragments_final.tsv.gz')
verbose = TRUE
cells<-names(x = cellmap)
Sys.time()
time_c<-Sys.time()
fragments <- Signac:::GetReadsInRegion(region = region, cellmap = cellmap,
cells = cells, tabix.file = tabix.file, verbose = verbose)########找出mapping到TSS position的片段
start.lookup <- start(x = region)
names(start.lookup) <- seq_along(region)
fragstarts <- start.lookup[fragments$ident] + 1

cut.df <- data.frame(position = c(fragments$start, fragments$end) - 
            fragstarts, cell = c(fragments$cell, fragments$cell), feature = c(fragments$ident, fragments$ident),
            stringsAsFactors = FALSE)######记录mapping片段到起始点和终止点的位置
cut.df <- cut.df[(cut.df$position > 0) & (cut.df$position <= 
            width(x = region)[[1]]), ]#########删除超覆盖的片段
tss_features<-1:max(fragments$ident)
Sys.time()
time_a<-Sys.time()
TSS_enrichment_frame_2K<-data.frame()
TSS_enrichment_frame_NDR<-data.frame()
for(i in 1:length(cells)){
temp_cut_df<-cut.df[cut.df$cell%in%cells[i],]
feature.vector <- seq_along(along.with = tss_features)
names(x = feature.vector) <- tss_features
feature.matrix.info <- feature.vector[temp_cut_df$feature]
feature_cut_matrix <- sparseMatrix(i = feature.matrix.info, j = temp_cut_df$position, 
            x = 1, dims = c(length(x = tss_features), width(x = region)[[1]]))

rownames(x = feature_cut_matrix) <- tss_features
colnames(x = feature_cut_matrix) <- seq_len(width(x = region)[[1]])
flanking.mean <- rowMeans(x = feature_cut_matrix[, c(2001:2101, 3901:4000)])
flanking.mean[is.na(x = flanking.mean)] <- 0
flanking.mean[flanking.mean == 0] <- mean(flanking.mean,na.rm = TRUE)
norm.matrix <- feature_cut_matrix/flanking.mean
TSS_enrichment_2K <- rowMeans(x = norm.matrix[, 2001:4000],na.rm = TRUE)
TSS_enrichment_NDR <- rowMeans(x = norm.matrix[, 2850:3050],na.rm = TRUE)
TSS_enrichment_frame_2K<-rbind(TSS_enrichment_frame_2K,t(TSS_enrichment_2K))
TSS_enrichment_frame_NDR<-rbind(TSS_enrichment_frame_NDR,t(TSS_enrichment_NDR))
print(round(i/length(cells),digits = 2))
}
TSS_enrichment_frame_2K<-t(TSS_enrichment_frame_2K)
colnames(TSS_enrichment_frame_2K)<-cells
TSS_enrichment_frame_NDR<-t(TSS_enrichment_frame_NDR)
colnames(TSS_enrichment_frame_NDR)<-cells
Sys.time()
time_b<-Sys.time()
print('for time')
time_b-time_a
print('all time')
time_b-time_c
save(TSS_enrichment_frame_2K,TSS_enrichment_frame_NDR,file='GSM5589392_pancreas_SM-IOBHS_rep1_TSS_enrichment.Rdata')

####################cell types of origin inferred by cell-free DNA  (CTORIC)
#######NMF去卷积,NNLS计算细胞比例
options(stringsAsFactors=F)
suppressMessages(require(Seurat))
suppressMessages(require(purrr))
suppressMessages(require(dplyr))
suppressMessages(require(tibble))
suppressMessages(require(NMF))
suppressMessages(require(Seurat))
suppressMessages(require(Matrix))
suppressMessages(require(dplyr))
suppressMessages(require(edgeR))
suppressMessages(require(SPOTlight))
set.seed(1234)
#####加载数据
load('/home/lijie/SPOT_ATAC/SPOT_test/tss_pro_table.Rdata')
load('/home/lijie/SPOT_ATAC/SPOT_test/cell_meta.Rdata')
ACINAR_index<-grep('ACINAR',cell_meta$cell_type)###删除ACINAR cell
cell_meta<-cell_meta[-ACINAR_index,]

ATAC_index1<-'ATAC1'######改变细胞编号
ATAC_cellmatrix1<-cell_meta[grep(ATAC_index1,cell_meta$cell_barcode),]######改变细胞编号
ATAC_index2<-'ATAC2'######改变细胞编号
ATAC_cellmatrix2<-cell_meta[grep(ATAC_index2,cell_meta$cell_barcode),]######改变细胞编号
ATAC_index3<-'ATAC3'######改变细胞编号
ATAC_cellmatrix3<-cell_meta[grep(ATAC_index3,cell_meta$cell_barcode),]######改变细胞编号
ATAC_index4<-'ATAC_n'######改变细胞编号
ATAC_cellmatrix4<-cell_meta[grep(ATAC_index4,cell_meta$cell_barcode),]######改变细胞编号
ATAC_index5<-'ATAC_p1'######改变细胞编号
ATAC_cellmatrix5<-cell_meta[grep(ATAC_index5,cell_meta$cell_barcode),]######改变细胞编号
ATAC_index6<-'ATAC_p2'######改变细胞编号
ATAC_cellmatrix6<-cell_meta[grep(ATAC_index6,cell_meta$cell_barcode),]######改变细胞编号
ATAC_index7<-'ATAC_p3'######改变细胞编号
ATAC_cellmatrix7<-cell_meta[grep(ATAC_index7,cell_meta$cell_barcode),]######改变细胞编号
ATAC_index8<-'ATAC_p4'######改变细胞编号
ATAC_cellmatrix8<-cell_meta[grep(ATAC_index8,cell_meta$cell_barcode),]######改变细胞编号
ATAC_cellmatrix<-rbind(ATAC_cellmatrix1,ATAC_cellmatrix2,ATAC_cellmatrix3,ATAC_cellmatrix4,ATAC_cellmatrix5,ATAC_cellmatrix6,
ATAC_cellmatrix7,ATAC_cellmatrix8)

load('/home/lijie/SPOT_ATAC/ATAC_1836591/ATAC_1836591_TSS_enrichment.Rdata')########选择scATAC TSS enrichment
new_TSS_enrichment_frame_NDR<-TSS_enrichment_frame_NDR
new_TSS_enrichment_frame_2K<-TSS_enrichment_frame_2K
load('/home/lijie/SPOT_ATAC/ATAC_T001835227/ATAC_T001835227_TSS_enrichment.Rdata')
new_TSS_enrichment_frame_NDR<-cbind(new_TSS_enrichment_frame_NDR,TSS_enrichment_frame_NDR)
new_TSS_enrichment_frame_2K<-cbind(new_TSS_enrichment_frame_2K,TSS_enrichment_frame_2K)
load('/home/lijie/SPOT_ATAC/ATAC_T001837519/ATAC_T001837519_TSS_enrichment.Rdata')
new_TSS_enrichment_frame_NDR<-cbind(new_TSS_enrichment_frame_NDR,TSS_enrichment_frame_NDR)
new_TSS_enrichment_frame_2K<-cbind(new_TSS_enrichment_frame_2K,TSS_enrichment_frame_2K)
load('/home/lijie/SPOT_ATAC/ATAC_T001837519_n/ATAC_T001837519_n_TSS_enrichment.Rdata')
new_TSS_enrichment_frame_NDR<-cbind(new_TSS_enrichment_frame_NDR,TSS_enrichment_frame_NDR)
new_TSS_enrichment_frame_2K<-cbind(new_TSS_enrichment_frame_2K,TSS_enrichment_frame_2K)
load('/home/lijie/SPOT_ATAC/GSM5589391_pancreas_SM-ADRUQ_rep1/GSM5589391_pancreas_SM-ADRUQ_rep1_TSS_enrichment.Rdata')
new_TSS_enrichment_frame_NDR<-cbind(new_TSS_enrichment_frame_NDR,TSS_enrichment_frame_NDR)
new_TSS_enrichment_frame_2K<-cbind(new_TSS_enrichment_frame_2K,TSS_enrichment_frame_2K)
load('/home/lijie/SPOT_ATAC/GSM5589392_pancreas_SM-IOBHS_rep1/GSM5589392_pancreas_SM-IOBHS_rep1_TSS_enrichment.Rdata')
new_TSS_enrichment_frame_NDR<-cbind(new_TSS_enrichment_frame_NDR,TSS_enrichment_frame_NDR)
new_TSS_enrichment_frame_2K<-cbind(new_TSS_enrichment_frame_2K,TSS_enrichment_frame_2K)
load('/home/lijie/SPOT_ATAC/GSM5589393_pancreas_SM-JF1NS_rep1/GSM5589393_pancreas_SM-JF1NS_rep1_TSS_enrichment.Rdata')
new_TSS_enrichment_frame_NDR<-cbind(new_TSS_enrichment_frame_NDR,TSS_enrichment_frame_NDR)
new_TSS_enrichment_frame_2K<-cbind(new_TSS_enrichment_frame_2K,TSS_enrichment_frame_2K)
load('/home/lijie/SPOT_ATAC/GSM5589394_pancreas_SM-JF1O6_rep1/GSM5589394_pancreas_SM-JF1O6_rep1_TSS_enrichment.Rdata')
new_TSS_enrichment_frame_NDR<-cbind(new_TSS_enrichment_frame_NDR,TSS_enrichment_frame_NDR)
new_TSS_enrichment_frame_2K<-cbind(new_TSS_enrichment_frame_2K,TSS_enrichment_frame_2K)
TSS_enrichment_frame_NDR<-new_TSS_enrichment_frame_NDR
rm(new_TSS_enrichment_frame_NDR)
TSS_enrichment_frame_2K<-new_TSS_enrichment_frame_2K
rm(new_TSS_enrichment_frame_2K)
rownames(TSS_enrichment_frame_2K)<-tss_pro_table$gene_name
rownames(TSS_enrichment_frame_NDR)<-tss_pro_table$gene_name


ATAC_cell<-unlist(lapply(strsplit(ATAC_cellmatrix$cell_barcode,'#'),function(x){return(x[[2]])}))
intersect_sample<-intersect(ATAC_cell,colnames(TSS_enrichment_frame_2K))
ATAC_cellmatrix<-ATAC_cellmatrix[match(intersect_sample,ATAC_cell),]
clust_vr=ATAC_cellmatrix$cell_type
TSS_enrichment_frame_NDR<-TSS_enrichment_frame_NDR[,intersect_sample]
TSS_enrichment_frame_2K<-TSS_enrichment_frame_2K[,intersect_sample]

setwd('/home/lijie/SPOT_ATAC/SPOT_test')
save(TSS_enrichment_frame_NDR,clust_vr,file='temp_scATAC_TSSscore_NDR.Rdata')
save(TSS_enrichment_frame_2K,clust_vr,file='temp_scATAC_TSSscore_2K.Rdata')


setwd('/home/lijie/SPOT_ATAC/SPOT_test')
load('temp_scATAC_TSSscore_NDR.Rdata')
load('temp_scATAC_TSSscore_2K.Rdata')
load('PDAC_sample_tss_150.Rdata')
colnames(temp_tss_150)[1:3]<-c('ATAC_1836591','ATAC_T001835227','ATAC_T001837519')
new_temp_tss_150<-temp_tss_150
load('healthy_sample_tss_150.Rdata')
new_temp_tss_150<-cbind(new_temp_tss_150,temp_tss_150)
temp_tss_150<-new_temp_tss_150
temp_tss_150[is.na(temp_tss_150)]<-0.01
temp_tss_150[temp_tss_150==0]<-0.01
tss_150_file<- -log2(temp_tss_150)
# tss_150_file<-tss_150_file-min(tss_150_file)#####让所有值都为正
load('PDAC_sample_tss_1000.Rdata')
colnames(temp_tss_1000)[1:3]<-c('ATAC_1836591','ATAC_T001835227','ATAC_T001837519')
new_temp_tss_1000<-temp_tss_1000
load('healthy_sample_tss_1000.Rdata')
new_temp_tss_1000<-cbind(new_temp_tss_1000,temp_tss_1000)
temp_tss_1000<-new_temp_tss_1000
temp_tss_1000[is.na(temp_tss_1000)]<-0.01
temp_tss_1000[temp_tss_1000==0]<-0.01
tss_1000_file<- -log2(temp_tss_1000)
# tss_1000_file<-tss_1000_file-min(tss_1000_file)#####让所有值都为正


#输入文件
load('/home/lijie/SPOT_ATAC/fragments/new/markerList_celltype.Rdata')######cutOff = "FDR <= 0.01 & Log2FC >= 2"
mtrx_spatial<-tss_1000_file
mtrx_spatial<-as.matrix(mtrx_spatial)
rownames(mtrx_spatial)<-tss_pro_table$gene_name

se_obj <- CreateSeuratObject(counts = TSS_enrichment_frame_2K)
se_obj@meta.data[["cell_type"]]<-clust_vr

######marker gene确定
# cl_n=100
# keep_ids <- lapply(split(se_obj@meta.data, se_obj@meta.data[["cell_type"]]), function(subdf) {
        # n_sample <- if_else(nrow(subdf) < cl_n, as.numeric(nrow(subdf)), 
            # as.numeric(cl_n))
		# set.seed(1234)
        # tmp_ds <- subdf[sample(seq_len(nrow(subdf)), n_sample), 
            # ] %>% tibble::rownames_to_column("barcodeID") %>% 
            # dplyr::pull(barcodeID)
        # return(tmp_ds)
    # }) %>% purrr::flatten_chr()
# se_obj<-se_obj[,keep_ids]
# select_frame<-TSS_enrichment_frame_2K[,keep_ids]
# cell_type_name<-names(table(se_obj@meta.data[["cell_type"]]))
# markerList<-list()
# for(i in 1:length(cell_type_name)){
# temp_markerList<-apply(select_frame,1,function(x){fc<-mean(x[colnames(select_frame)[se_obj@meta.data[["cell_type"]]==cell_type_name[i]]])/mean(x[colnames(select_frame)[!se_obj@meta.data[["cell_type"]]==cell_type_name[i]]]);
# p_value<-wilcox.test(x[colnames(select_frame)[se_obj@meta.data[["cell_type"]]==cell_type_name[i]]],x[colnames(select_frame)[!se_obj@meta.data[["cell_type"]]==cell_type_name[i]]]);return(c(fc,p_value$p.value))})
# markerList[[i]]<-temp_markerList
# print(i)
# }
# names(markerList)<-cell_type_name
# #cut_off
# cutoff_fc<-2
# cutoff_p<-0.05
# cell_type_marker<-data.frame()
# for(i in 1:length(markerList)){
# temp_cell_type_marker<-t(markerList[[i]])
# temp_frame<-temp_cell_type_marker[temp_cell_type_marker[,1]>=cutoff_fc&temp_cell_type_marker[,2]<=cutoff_p,]
# temp_frame<-as.data.frame(na.omit(temp_frame))
# temp_complete_frame<-data.frame(name=rownames(temp_frame),temp_frame,cell_type=rep(names(markerList)[i],dim(temp_frame)[1]))
# colnames(temp_complete_frame)<-c('gene','avg_logFC','p_val','cluster')
# cell_type_marker<-rbind(cell_type_marker,temp_complete_frame)
# }
# save(cell_type_marker,markerList,file='TSS_score_cell_type_marker.Rdata')

# load('TSS_score_cell_type_marker.Rdata')
#######
cell_type_marker<-rbind(data.frame(markerList$B[,c('name','Log2FC','FDR')],cluster=rep('B',dim(markerList$B)[1])),
data.frame(markerList$DC[,c('name','Log2FC','FDR')],cluster=rep('DC',dim(markerList$DC)[1])),
data.frame(markerList$DUC[,c('name','Log2FC','FDR')],cluster=rep('DUC',dim(markerList$DUC)[1])),
data.frame(markerList$DUC_n[,c('name','Log2FC','FDR')],cluster=rep('DUC_n',dim(markerList$DUC_n)[1])),
data.frame(markerList$ENDOCRINE[,c('name','Log2FC','FDR')],cluster=rep('ENDOCRINE',dim(markerList$ENDOCRINE)[1])),
data.frame(markerList$ENDOTHELIAL[,c('name','Log2FC','FDR')],cluster=rep('ENDOTHELIAL',dim(markerList$ENDOTHELIAL)[1])),
data.frame(markerList$FIBROBLAST[,c('name','Log2FC','FDR')],cluster=rep('FIBROBLAST',dim(markerList$FIBROBLAST)[1])),
data.frame(markerList$MAC[,c('name','Log2FC','FDR')],cluster=rep('MAC',dim(markerList$MAC)[1])),
data.frame(markerList$NKT[,c('name','Log2FC','FDR')],cluster=rep('NKT',dim(markerList$NKT)[1])),
data.frame(markerList$T[,c('name','Log2FC','FDR')],cluster=rep('T',dim(markerList$T)[1])))
colnames(cell_type_marker)<-c('gene','avg_logFC','p_val','cluster')

# ######寻找DUC marker,和T,B等细胞marker
# rownames(tss_1000_file)<-tss_pro_table$gene_name
# tss_1000_file<-as.data.frame(tss_1000_file)
# #DUC
# temp_frame<-na.omit(tss_1000_file[cell_type_marker[cell_type_marker$cluster%in%'DUC',]$gene,])
# DUC_pvalue<-apply(temp_frame,1,function(x){a<-wilcox.test(x[1:39],x[40:139],alternative='great');return(a$p.value)})
# DUC_marker<-names(DUC_pvalue[which(DUC_pvalue<0.1)])####卡的阈值
# cell_type_DUC<-cell_type_marker[cell_type_marker$gene%in%DUC_marker,]
# cell_type_DUC<-cell_type_DUC[cell_type_DUC$cluster=='DUC',]

# #T
# extra_cell_type<-c('B','DC','DUC_n','ENDOCRINE','ENDOTHELIAL','FIBROBLAST','MAC','NKT','T')
# new_cell_type_marker<-data.frame()
# for(i in 1:length(extra_cell_type)){
# temp_frame<-na.omit(tss_1000_file[cell_type_marker[cell_type_marker$cluster%in%extra_cell_type[i],]$gene,])
# T_pvalue<-apply(temp_frame,1,function(x){a<-wilcox.test(x[40:90],x[91:139]);return(a$p.value)})
# T_marker<-names(T_pvalue[which(T_pvalue>0.05)])####卡的阈值
# cell_type_T<-cell_type_marker[cell_type_marker$gene%in%T_marker,]
# cell_type_T<-cell_type_T[cell_type_T$cluster==extra_cell_type[i],]
# new_cell_type_marker<-rbind(new_cell_type_marker,cell_type_T)
# }
# cell_type_marker<-rbind(cell_type_DUC,new_cell_type_marker)

# PDAC_W<-tss_1000_file[cell_type_DUC$gene,1:39]
# healthy_W<-tss_1000_file[cell_type_DUC$gene,40:139]

#########downsample
###########marker gene and hvg确定,downsample cl_n=100 clust_vr=cell_type
method = "nsNMF"
hvg=3000
cl_n=100

# se_obj <- Seurat::FindVariableFeatures(object = se_obj, 
            # nfeatures = hvg)
# keep_genes <- unique(c(VariableFeatures(se_obj), cell_type_marker$gene))
keep_genes <- unique(cell_type_marker$gene)

keep_ids <- lapply(split(se_obj@meta.data, se_obj@meta.data[["cell_type"]]), function(subdf) {
        n_sample <- if_else(nrow(subdf) < cl_n, as.numeric(nrow(subdf)), 
            as.numeric(cl_n))
		set.seed(164)#j=164
        tmp_ds <- subdf[sample(seq_len(nrow(subdf)), n_sample), 
            ] %>% tibble::rownames_to_column("barcodeID") %>% 
            dplyr::pull(barcodeID)
        return(tmp_ds)
    }) %>% purrr::flatten_chr()
se_obj <- se_obj[keep_genes, keep_ids]

######train nmf
slot = "counts"
assay = "RNA"
transf = "uv"
cluster_markers<-cell_type_marker
########初始NMF W与H矩阵设定
ntop= NULL
seed_init_mtrx_nmf<-function (cluster_markers, se_obj, ntop = NULL) 
{
    se_nmf_ready <- prep_seobj_topic_fun(se_obj = se_obj)
    k <- length(unique(se_obj@meta.data[["cell_type"]]))
    if (is.null(ntop)) 
        ntop <- max(table(cluster_markers$cluster))
    cluster_markers_cut <- suppressMessages(cut_markers2(markers = cluster_markers, 
        ntop = ntop))
    cluster_markers_uniq <- lapply(unique(cluster_markers_cut$cluster), 
        function(clust) {
            ls1 <- cluster_markers_cut[cluster_markers_cut$cluster == 
                clust, "gene"]
            ls2 <- cluster_markers_cut[cluster_markers_cut$cluster != 
                clust, "gene"]
            ls1_unique <- ls1[!ls1 %in% ls2]
            return(cluster_markers_cut[cluster_markers_cut$cluster == 
                clust & cluster_markers_cut$gene %in% ls1_unique, 
                ])
        }) %>% bind_rows()#####取每类里的唯一marker
    seedgenes <- matrix(nrow = k, ncol = ncol(se_nmf_ready), 
        data = 1e-10)
    colnames(seedgenes) <- colnames(se_nmf_ready)
    for (i in seq_len(k)) {
        clust_row <- cluster_markers_uniq$cluster == as.character(unique(se_obj@meta.data[["cell_type"]])[[i]])
        seedgenes[i, as.character(cluster_markers_uniq[clust_row, 
            "gene"])] = cluster_markers_uniq[clust_row, "weight"]
    }###########cell type marker gene权重设置为weight = 1 - p_val，其他为1e-10
    W <- t(seedgenes)
    H <- matrix(data = 1e-10, nrow = k, ncol = nrow(se_nmf_ready))
    for (i in seq_len(nrow(se_nmf_ready))) {
        h_row <- which(unique(se_obj@meta.data[["cell_type"]]) == 
            se_obj@meta.data[["cell_type"]][i])
        H[h_row, i] <- 1
    }###########用1标记出哪些是具体的细胞，1e-10代表不是具体的那个细胞
    rownames(W) <- rownames(se_obj@assays$RNA@counts)
    colnames(H) <- colnames(se_obj@assays$RNA@counts)
    return(list(W = W, H = H))
}
###
print("Preparing Gene set")
mtrx_sc <- as.matrix(Seurat::GetAssayData(se_obj, assay = assay, slot = slot))######刚下采样的矩阵
genes_0_sc <- which(!rowSums(mtrx_sc == 0) == ncol(mtrx_sc))#####删除整行为0
cell_0_sc <- which(!colSums(mtrx_sc == 0) == nrow(mtrx_sc))#####删除整列为0
se_obj <- se_obj[genes_0_sc,cell_0_sc]
genes_0_sp <- which(!rowSums(as.matrix(mtrx_spatial) == 0) == 
        ncol(mtrx_spatial))
mtrx_spatial <- mtrx_spatial[genes_0_sp, ]
genes_spatial <- rownames(mtrx_spatial)
genes_sc <- rownames(Seurat::GetAssayData(se_obj, assay = assay, 
        slot = slot))
if (length(intersect(genes_sc, genes_spatial)) < 10) 
        stop("Not enough genes in common between the single-cell and mixture dataset.")
se_obj <- se_obj[intersect(genes_sc, genes_spatial), ]
mtrx_sc <- as.matrix(Seurat::GetAssayData(se_obj, assay = assay, 
        slot = slot))
cluster_markers <- cluster_markers[cluster_markers$gene%in%rownames(se_obj), ]
print("Normalizing count matrix")
if (transf == "uv") {
        count_mtrx_t <- scale(t(mtrx_sc), center = FALSE, scale = apply(mtrx_sc, 
            1, sd, na.rm = TRUE))###仅进行标准差标准化
        count_mtrx <- t(count_mtrx_t)
    }else if (transf == "raw") {
        count_mtrx <- mtrx_sc
    }
cell_0_sc <- which(!colSums(count_mtrx == 0) == nrow(count_mtrx))#####删除整列为0
se_obj <- se_obj[,cell_0_sc]
count_mtrx<-count_mtrx[,cell_0_sc]
k <- length(unique(se_obj@meta.data[["cell_type"]]))
start_t <- Sys.time()
if (method == "nsNMF") 
        mod <- "NMFns"
if (is.numeric(hvg)) {
        print("Seeding initial matrices")
        init_mtrx <- seed_init_mtrx_nmf(cluster_markers = cluster_markers, 
            se_obj = se_obj, ntop = NULL)#####采用下采样矩阵进行初始矩阵W,H分解
        nmf_init <- NMF::nmfModel(W = init_mtrx[["W"]], H = init_mtrx[["H"]], 
            model = mod)
        print("Training...")
        nmf_mod <- NMF::nmf(x = count_mtrx, rank = k, seed = nmf_init, 
            method = method)
}
total_t <- round(difftime(Sys.time(), start_t, units = "mins"), 
        2)
print(sprintf("Time to train NMF model was %smins", total_t))
nmf_mod_ls<-list(nmf_mod, as.vector(se_obj@meta.data[["cell_type"]]))


###############返回每种细胞类型的H中位数矩阵
h = coef(nmf_mod_ls[[1]])
train_cell_clust = nmf_mod_ls[[2]]
suppressMessages(require(tibble))
suppressMessages(require(dplyr))
h_ds <- data.frame(t(h))
h_ds[, "clust_vr"] <- train_cell_clust
ct_topic_profiles <- h_ds %>% dplyr::group_by(clust_vr) %>% 
dplyr::summarise_all(list(median)) %>% tibble::remove_rownames() %>% 
tibble::column_to_rownames(var = "clust_vr") %>% as.matrix()######对每个细胞类型取中值
ct_topic_profiles_t <- t(ct_topic_profiles)
colnames(ct_topic_profiles_t) <- gsub("[[:punct:]]|[[:blank:]]", 
".", colnames(ct_topic_profiles_t))
ct_topic_profiles<-ct_topic_profiles_t

save(nmf_mod_ls,ct_topic_profiles,file='TSS_2K_ATACallpublished_processed_W.Rdata')
########NNLS
mixture_transcriptome = mtrx_spatial
nmf_mod = nmf_mod_ls[[1]]
transf = "uv"
reference_profiles = ct_topic_profiles
min_cont = 0.01
#######以W为基础，获取空转的系数矩阵
predict_spatial_mixtures_nmf<-function (nmf_mod, mixture_transcriptome, transf) 
{
    suppressMessages(require(nnls))
    suppressMessages(require(edgeR))
    keep_genes <- rownames(basis(nmf_mod))[rownames(basis(nmf_mod)) %in% 
        rownames(mixture_transcriptome)]
    mixture_transcriptome_subs <- as.matrix(mixture_transcriptome[keep_genes, 
        ])
	#index_var<-which(apply(mixture_transcriptome_subs, 1, sd,na.rm = TRUE)==0)
    # if (transf == "uv") {
        # count_mtrx <- scale(t(mixture_transcriptome_subs), center = FALSE, 
            # scale = apply(mixture_transcriptome_subs, 1, sd, 
                # na.rm = TRUE))###仅进行标准差标准化
        # count_mtrx <- t(count_mtrx)
        # pos_0 <- which(rowSums(is.na(count_mtrx)) == ncol(count_mtrx))
        # count_mtrx[pos_0, ] <- 0
    # }
    #count_mtrx<-count_mtrx[-index_var,]
    count_mtrx<-mixture_transcriptome_subs
	W <- basis(nmf_mod)
	#W<-W[-index_var,]
    coef_pred <- matrix(data = NA, nrow = ncol(W), ncol = ncol(count_mtrx))
    colnames(coef_pred) <- colnames(count_mtrx)
    for (i in seq_len(ncol(count_mtrx))) {
        nnls_pred <- nnls::nnls(A = W, b = count_mtrx[, i])
        coef_pred[, i] <- nnls_pred$x
    }
    return(coef_pred)
}
# mixture_deconvolution_nmf######以NNLS[nnls(A=x,B=x)]为基础，计算出每种细胞的比例

    suppressMessages(require(nnls))
    profile_mtrx <- predict_spatial_mixtures_nmf(nmf_mod = nmf_mod, 
        mixture_transcriptome = mixture_transcriptome, transf = transf)
	profile_mtrx[3,40:139]<-0
	profile_mtrx <-coef_pred
	# profile_mtrx[3,]<-profile_mtrx[3,]*100^new_W
    decon_mtrx <- matrix(data = NA, nrow = ncol(profile_mtrx), 
        ncol = ncol(reference_profiles) + 1)
    colnames(decon_mtrx) <- c(colnames(reference_profiles), "res_ss")
    print("Deconvoluting spots")
    total <- ncol(profile_mtrx)
    pb <- txtProgressBar(min = 0, max = total, style = 3)
    for (i in seq_len(ncol(profile_mtrx))) {
        nnls_pred <- nnls::nnls(A = reference_profiles, b = profile_mtrx[, 
            i])
        weights <- nnls_pred$x
        comp <- weights/sum(weights)
        comp[comp < min_cont] <- 0
        weights[comp < min_cont] <- 0
        comp_prop <- comp/sum(comp)
        comp_prop[is.na(comp_prop)] <- 0
        fit_null <- 0
        tot_ss <- sum((profile_mtrx[, i] - fit_null)^2)
        unexpl_ss <- nnls_pred$deviance/tot_ss
        decon_mtrx[i, 1:(ncol(decon_mtrx) - 1)] <- comp_prop
        decon_mtrx[i, ncol(decon_mtrx)] <- unexpl_ss
        setTxtProgressBar(pb, i)
    }
    close(pb)
	rownames(decon_mtrx)<-colnames(mtrx_spatial)
    View(decon_mtrx)
save(decon_mtrx,file='decon_mtrx.Rdata')
length(which(decon_mtrx[40:139,3]==0))

####细胞比例barplot
temp_cell_fraction<-c()
for(i in 1:10){
temp_cell_fraction<-c(temp_cell_fraction,decon_mtrx[,i])
}

df_fc<-data.frame(cell_fraction=temp_cell_fraction,cell_type=factor(c(rep(colnames(decon_mtrx)[1:10],each=139)),levels=colnames(decon_mtrx)[1:10]),
type=factor(c(rep(c(rep("tumor",39),rep("normal",100)),10)),levels=c('tumor','normal')))

df_fc$combine_name<-factor(paste(df_fc$gene,df_fc$type),levels=c("CFTR tumor","CFTR normal","ZNF616 tumor","ZNF616 normal","NFIA tumor",
"NFIA normal","MDFI tumor","MDFI normal","ZNF667 tumor","ZNF667 normal",
"MYO10 tumor","MYO10 normal","ZFPL1 tumor","ZFPL1 normal","ENSA tumor","ENSA normal","PIFO tumor","PIFO normal","GGA1 tumor","GGA1 normal"))



library(ggplot2)
library(ggpubr)
color = c("#BE4E4D","#35A1D3")
pdf('all_cell_fraction.pdf',width=20,height=10)
my_comparisons <- list(c("tumor", "normal"))
df_fc %>% ggplot(aes(x=cell_type,y=cell_fraction,fill=type))+
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
  stat_compare_means(label.y = 0.5)+ylim(0,0.8)
dev.off()


#########细胞比例的pie plot
options(stringsAsFactors=F)
library(ggplot2)
library(dplyr)

# Create Data
count.matrix<-rbind(t(colnames(new_mtr)),new_mtr)

count.data <- data.frame(class = c('STAD','Healthy','Gastritis','BRCA','COAD','PDAC','ESCA','LUAD','THCA','ACCx','OV','UCEC','LIHC','GBM',
'Bone disease','Cerebral hemorrhage','Heart disease','Benigh tumor,cyst','Lithiasis','Lleus','Other'),n = c(120,145,83,11,13,40,1,3,3,2,1,2,1,1,14,
3,6,17,8,2,8))
count.data$prop<-signif(count.data$n/sum(count.data$n),digits = 2)*100
# Compute the position of labels

count.data<-count.matrix[,c(1,4)]
count.data<-as.data.frame(count.data)
count.data[,2]<-as.numeric(count.data[,2])
colnames(count.data)[2]<-'prop'

count.data <- count.data %>%
arrange(desc(class)) %>%
mutate(lab.ypos = cumsum(prop) - 0.5*prop)

count.data$prop<-signif(count.data$prop,digits = 2)

# Donut chart
ggplot(count.data, aes(x = 2, y = prop, fill = class)) +
  geom_bar(stat = "identity", color = "white") +
  coord_polar(theta = "y", start = 0)+
  theme_void()+
  geom_text(aes(y = lab.ypos, label = prop), color = "white", size=6) +
  xlim(0.5, 2.5)


#########细胞比例的pie plot
options(stringsAsFactors=F)
library(ggplot2)
library(dplyr)

# Create Data
i<-3
count.data <- data.frame(class = colnames(decon_mtrx)[1:10],n = decon_mtrx[i,1:10])
count.data$prop<-signif(count.data$n,digits = 2)*100
# Compute the position of labels

count.data<-count.data[,c(1,3)]
count.data<-as.data.frame(count.data)
count.data[,2]<-as.numeric(count.data[,2])
colnames(count.data)[2]<-'prop'

count.data <- count.data %>%
arrange(desc(class)) %>%
mutate(lab.ypos = cumsum(prop) - 0.5*prop)

count.data$prop<-signif(count.data$prop,digits = 2)

# Donut chart
ggplot(count.data, aes(x = 2, y = prop, fill = class)) +
  geom_bar(stat = "identity", color = "white") +
  coord_polar(theta = "y", start = 0)+
  theme_void()+
  geom_text(aes(y = lab.ypos, label = prop), color = "white", size=6) +
  xlim(0.5, 2.5)

#######cfDNA infer与sc PBMC的相关性
a<-apply(decon_mtrx,2,median)
cfDNA_infer_value<-c(a[c(10,9,8,1,2)])
pb_5931<-c(1263,633,276,312,11)
pb_6207<-c(1203,668,657,43,10)
cor.test(pb_5931,c(a[c(10,9,8,1,2)]))
cor.test(pb_6207,c(a[c(10,9,8,1,2)]))
dfm<-data.frame(sc_pbmc=c(pb_5931,pb_6207),cfDNA_infer=c(cfDNA_infer_value,cfDNA_infer_value),group=c(rep('pb_5931',length(pb_5931)),
rep('pb_6207',length(pb_6207))))

library(ggplot2)
pdf('pbmc_cfDNA_infer_value_cor.pdf',width=11,height=10)
p<-ggplot(dfm,aes(x=sc_pbmc,y=cfDNA_infer,color=group,shape=group))+geom_smooth(method='lm')+geom_point()+ylim(0,max(0.4))
p+theme_bw()+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())
dev.off()

dfm<-data.frame(sc_pbmc=pb_6207,cfDNA_infer=cfDNA_infer_value,group=c(rep('sc_pbmc',length(pb_6207)),
rep('cfDNA_infer',length(cfDNA_infer_value))))

library(ggplot2)
pdf('pb_6207_cfDNA_infer_value_cor.pdf',width=11,height=10)
p<-ggplot(dfm,aes(x=sc_pbmc,y=cfDNA_infer,fill=group))+geom_smooth(method='lm')+geom_point(color='#F0A13E')
p+theme_bw()+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())
dev.off()
