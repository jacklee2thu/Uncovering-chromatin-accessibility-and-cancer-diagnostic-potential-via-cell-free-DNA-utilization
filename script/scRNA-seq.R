########gastric single cell为细胞培养做准备
options(stringsAsFactors=F)
setwd('/Share2/home/lanxun3/jacklee/cfDNA/filed_stanford_gastic_singlecell/count')
sample_all<-dir()
tumor_sample<-sample_all[grep('t',sample_all)]
normal_sample<-sample_all[grep('n',sample_all)]
pbmc_sample<-sample_all[grep('pbmc',sample_all)]

temp_mat<-data.frame()
temp_sample<-pbmc_sample
library(Matrix)
for(i in 1:length(temp_sample)){
matrix_dir = paste0(temp_sample[i],'/')
barcode.path <- paste0(matrix_dir, "barcodes.tsv.gz")
features.path <- paste0(matrix_dir, "features.tsv.gz")
matrix.path <- paste0(matrix_dir, "matrix.mtx.gz")
mat <- readMM(file = matrix.path)
feature.names = read.delim(features.path, 
                           header = FALSE,
                           stringsAsFactors = FALSE)
barcode.names = read.delim(barcode.path, 
                           header = FALSE,
                           stringsAsFactors = FALSE)
colnames(mat) = barcode.names$V1
rownames(mat) = feature.names$V2
c_mat<-as.matrix(mat)
c_mat<-as.data.frame(c_mat)
colnames(c_mat)<-paste0(temp_sample[i],'_',colnames(c_mat))

temp_mat<-rbind(temp_mat,t(c_mat))
}
tumor_mat<-t(temp_mat)
normal_mat<-t(temp_mat)
pbmc_mat<-t(temp_mat)

save(pbmc_mat,file='pbmc_mat.Rdata')


#####细胞分类
options(stringsAsFactors=F)
setwd("/Share2/home/lanxun3/jacklee/cfDNA/filed_stanford_gastic_singlecell")
load('tumor_mat.Rdata')
# load('pbmc_mat.Rdata')
load('normal_mat.Rdata')
setwd('/Share2/home/lanxun3/jacklee/cfDNA/filed_stanford_gastic_singlecell/data_precess')
all_sample_file<-cbind(tumor_mat,normal_mat)
library(dplyr)
library(Seurat)
all_sample <- CreateSeuratObject(counts = all_sample_file, project = "all_sample",min.cells = 3,min.features = 200)
all_sample[["percent.mt"]] <- PercentageFeatureSet(all_sample, pattern = "^MT-")
head(all_sample@meta.data, 5)
all_sample <- subset(all_sample, subset = nFeature_RNA > 200 & percent.mt < 20 )#####此处没有限制4000上限

dim(all_sample@assays$RNA)
#save filtered object

save(all_sample,file = "all_sample_qc.Rdata")

######细胞分类,findmarker
options(stringsAsFactors=F)
setwd("/Share2/home/lanxun3/jacklee/cfDNA/filed_stanford_gastic_singlecell/data_precess")
load('all_sample_qc.Rdata')
library(Seurat)
all_sample_normalized <- NormalizeData(all_sample,verbose = F)
all_sample_feature <- FindVariableFeatures(all_sample_normalized,selection.method = 'vst',nfeatures = 2000,verbose = F)
#######################################################
feature_var <- all_sample_feature@assays$RNA@var.features
cell_cycle_gene<-c("ASPM","CENPE","CENPF","DLGAP5","MKI67","NUSAP1","PCLAF","STMN1","TOP2A","TUBB")
feature_var <- setdiff(feature_var,cell_cycle_gene)

#######################################################
all_sample_scale <- ScaleData(all_sample_feature,verbose = F)
all_sample_pca <- RunPCA(all_sample_scale,features = feature_var,verbose = F)
pdf('cluster_num.pdf')
ElbowPlot(all_sample_pca)
dev.off()
set.seed(12345678)
all_sample_pca <- FindNeighbors(all_sample_pca,dims = 1:15)
all_sample_pca <- FindClusters(all_sample_pca,resolution = 0.8)
sample_cluster<-Idents(all_sample_pca)#####样本和cluster对应
sample_summary<-data.frame(names(sample_cluster),sample_cluster)
colnames(sample_summary)<-c('sample_name','cluster')
table(sample_summary$cluster)
save(sample_summary,file='sample_summary.Rdata')


all_sample_umap <- RunUMAP(all_sample_pca,dims = 1:15)
pdf('gastric_cancer_UMAP_scRNA_cluster.pdf')
DimPlot(object = all_sample_umap, reduction = "umap",label = T)+ NoLegend()
dev.off()
all_sample_umap@meta.data$sample_name<-strtrim(sample_summary$sample_name,7)
all_sample_umap@meta.data$sample_type<-strtrim(lapply(strsplit(sample_summary$sample_name,'_'),function(x){x[[2]]}),1)

pdf('gastric_cancer_UMAP_scRNA_sample.pdf',width=11,height=10)
DimPlot(object = all_sample_umap, reduction = "umap",group.by = 'sample_name')
dev.off()
pdf('gastric_cancer_UMAP_scRNA_sample_type.pdf')
DimPlot(object = all_sample_umap, reduction = "umap",group.by = 'sample_type',cols=c('n'='#447CAA','t'='#841B1F'))+ NoLegend()
dev.off()
save(all_sample_umap,file="all_sample_umap.Rdata")
#############################findmarker
options(stringsAsFactors=F)
setwd("/Share2/home/lanxun3/jacklee/cfDNA/filed_stanford_gastic_singlecell/data_precess")
load('all_sample_umap.Rdata')
library(Seurat)
library(dplyr)
#cluster_0.markers <- FindMarkers(object = all_sample_umap, ident.1 = 0, min.pct = 0.25,only.pos = T)
#pdac_M0 <- head(cluster_0.markers, n = 5)

##marker gene heatmap
c_spatial.markers <- FindAllMarkers(all_sample_umap, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
save(c_spatial.markers,file='c_spatial_markers.Rdata')
load('c_spatial_markers.Rdata')
top10 <- c_spatial.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
pdf('marker_gene.pdf')
DoHeatmap(all_sample_umap, features = top10$gene) + NoLegend()
dev.off()

########可视化特征基因
load('c_spatial_markers.Rdata')
feature_gene<-c('KRT18','EPCAM','CD3D','CD68','CHGA','CDH5','DCN','ID2','MS4A1','CD79A','IGLL5','MS4A2')
pdf('feature_gene_MEIS2.pdf')
VlnPlot(all_sample_umap,'MEIS2',pt.size = -1)
dev.off()
pdf('feature_gene_2.pdf',width=10,height=12)
VlnPlot(all_sample_umap,features=feature_gene,pt.size = -1)
dev.off()


c_spatial.markers[c_spatial.markers$gene=='CDH5',]$cluster

new.cluster.ids <- c('T','T','epithelium','fibroblast','epithelium','epithelium','B','B','epithelium','epithelium','endothelial','epithelium',
'macrophage','epithelium','epithelium','T','fibroblast','epithelium','macrophage','epithelium','T','fibroblast','T','endothelial','mast','macrophage',
'epithelium','endocrine','plasma','endocrine','epithelium','B','endothelial','epithelium','macrophage')
					 
names(new.cluster.ids) <- levels(all_sample_umap)

all_sample_umap <- RenameIdents(all_sample_umap, new.cluster.ids)
pdf('gastric_cancer_UMAP_scRNA_cell_types.pdf')
DimPlot(object = all_sample_umap, reduction = "umap",label = F,label.size = 4,repel=T,cols=c('B'='#1AB4B8','T'='#ED877F','endocrine'='#E06AA4',
'plasma'='#008B00','mast'='#8585B9','epithelium'='#DDB056','macrophage'="#7DC05B",'fibroblast'="#7DC6AC",'endothelial'="#8B3A3A"))+ NoLegend()
dev.off()

sample_cell<-Idents(all_sample_umap)
cell_summary<-data.frame(strtrim(names(sample_cell),7),sample_cell)
colnames(cell_summary)<-c('sample_name','cell_type')
save(cell_summary,file='cell_summary.Rdata')

pdf('feature_gene_dot.pdf',width=10,height=10)
all_sample_umap@active.ident<-factor(all_sample_umap@active.ident,levels=rev(c('epithelium','T','macrophage',
'endocrine','endothelial','fibroblast','DC','B','plasma','mast')))
DotPlot(all_sample_umap,features=rev(feature_gene),cols=c('#BEBEBE','#FF0000'))
dev.off()
########marker gene Dot plot
options(stringsAsFactors=F)
setwd("/Share2/home/lanxun3/jacklee/cfDNA/filed_stanford_gastic_singlecell/data_precess")
load('all_sample_umap.Rdata')
load('cell_summary.Rdata')
epithelial_cell<-grep('epithelium',cell_summary$cell_type)
library(Seurat)
library(dplyr)
markerGenes<-c('RAE1','GPR85','CALU','AK9','TRAF2','TENT2','ABCA10','LRRIQ1','SPECC1','MEIS2')
new_sample_umap<-subset(all_sample_umap[,epithelial_cell])
new_sample_umap@active.ident<-factor(new_sample_umap@meta.data$sample_type,levels=c('n','t'))

pdf('cfDNA_marker_gene_dot.pdf',width=10,height=5)
DotPlot(new_sample_umap,features=rev(markerGenes),cols=c('#BEBEBE','#FF0000'),scale=F)
dev.off()

pdf('feature_gene_TENT2.pdf')
FeaturePlot(all_sample_umap,features='TENT2',cols=c('#BEBEBE','#FF0000'),order=T,max.cutoff=0.1)+ NoLegend()
dev.off()



#####PBMC细胞分类,比例
options(stringsAsFactors=F)
setwd("/Share2/home/lanxun3/jacklee/cfDNA/filed_stanford_gastic_singlecell")
load('pbmc_mat.Rdata')
setwd('/Share2/home/lanxun3/jacklee/cfDNA/filed_stanford_gastic_singlecell/pbmc')
pbmc_sample_file<-pbmc_mat
library(dplyr)
library(Seurat)
pbmc_sample <- CreateSeuratObject(counts = pbmc_sample_file, project = "pbmc_sample",min.cells = 3,min.features = 200)
pbmc_sample[["percent.mt"]] <- PercentageFeatureSet(pbmc_sample, pattern = "^MT-")
head(pbmc_sample@meta.data, 5)
pbmc_sample <- subset(pbmc_sample, subset = nFeature_RNA > 200 & percent.mt < 20 )#####此处没有限制4000上限

dim(pbmc_sample@assays$RNA)
#save filtered object

save(pbmc_sample,file = "pbmc_sample_qc.Rdata")

######细胞分类,findmarker
options(stringsAsFactors=F)
setwd("/Share2/home/lanxun3/jacklee/cfDNA/filed_stanford_gastic_singlecell/pbmc")
load('pbmc_sample_qc.Rdata')
library(Seurat)
pbmc_sample_normalized <- NormalizeData(pbmc_sample,verbose = F)
pbmc_sample_feature <- FindVariableFeatures(pbmc_sample_normalized,selection.method = 'vst',nfeatures = 2000,verbose = F)
#######################################################
feature_var <- pbmc_sample_feature@assays$RNA@var.features
cell_cycle_gene<-c("ASPM","CENPE","CENPF","DLGAP5","MKI67","NUSAP1","PCLAF","STMN1","TOP2A","TUBB")
feature_var <- setdiff(feature_var,cell_cycle_gene)

#######################################################
pbmc_sample_scale <- ScaleData(pbmc_sample_feature,verbose = F)
pbmc_sample_pca <- RunPCA(pbmc_sample_scale,features = feature_var,verbose = F)

set.seed(12345678)
pbmc_sample_pca <- FindNeighbors(pbmc_sample_pca,dims = 1:15)
pbmc_sample_pca <- FindClusters(pbmc_sample_pca,resolution = 0.8)
sample_cluster<-Idents(pbmc_sample_pca)#####样本和cluster对应
sample_summary<-data.frame(names(sample_cluster),sample_cluster)
colnames(sample_summary)<-c('sample_name','cluster')
table(sample_summary$cluster)
save(sample_summary,file='sample_summary.Rdata')


pbmc_sample_umap <- RunUMAP(pbmc_sample_pca,dims = 1:15)
pdf('pbmc_UMAP_scRNA_cluster.pdf')
DimPlot(object = pbmc_sample_umap, reduction = "umap",label = T)+ NoLegend()
dev.off()
pbmc_sample_umap@meta.data$sample_name<-strtrim(sample_summary$sample_name,7)
pbmc_sample_umap@meta.data$sample_type<-strtrim(lapply(strsplit(sample_summary$sample_name,'_'),function(x){x[[2]]}),1)

pdf('pbmc_UMAP_scRNA_sample.pdf',width=11,height=10)
DimPlot(object = pbmc_sample_umap, reduction = "umap",group.by = 'sample_name')
dev.off()
pdf('pbmc_UMAP_scRNA_sample_type.pdf')
DimPlot(object = pbmc_sample_umap, reduction = "umap",group.by = 'sample_type',cols=c('n'='#447CAA','t'='#841B1F'))+ NoLegend()
dev.off()
save(pbmc_sample_umap,file="pbmc_sample_umap.Rdata")
#############################findmarker
options(stringsAsFactors=F)
setwd("/Share2/home/lanxun3/jacklee/cfDNA/filed_stanford_gastic_singlecell/pbmc")
load('pbmc_sample_umap.Rdata')
library(Seurat)
library(dplyr)
#cluster_0.markers <- FindMarkers(object = all_sample_umap, ident.1 = 0, min.pct = 0.25,only.pos = T)
#pdac_M0 <- head(cluster_0.markers, n = 5)

##marker gene heatmap
c_spatial.markers <- FindAllMarkers(pbmc_sample_umap, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
save(c_spatial.markers,file='c_spatial_markers.Rdata')
load('c_spatial_markers.Rdata')
top10 <- c_spatial.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
pdf('marker_gene.pdf')
DoHeatmap(pbmc_sample_umap, features = top10$gene) + NoLegend()
dev.off()

########可视化特征基因
load('c_spatial_markers.Rdata')
feature_gene<-c('CD3D','MS4A1','CD79A','NKG7','CD68','AIF1','FCER1A')
pdf('feature_gene_MEIS2.pdf')
VlnPlot(all_sample_umap,'MEIS2',pt.size = -1)
dev.off()
pdf('feature_gene_2.pdf',width=10,height=12)
VlnPlot(all_sample_umap,features=feature_gene,pt.size = -1)
dev.off()


c_spatial.markers[c_spatial.markers$gene=='CD3D',]$cluster

new.cluster.ids <- c('T','NK','NK','mac','T','T','T','mac','B','T','T','T','T','B','mac','NK','mac','T','DC')
					 
names(new.cluster.ids) <- levels(pbmc_sample_umap)

pbmc_sample_umap <- RenameIdents(pbmc_sample_umap, new.cluster.ids)
pdf('pbmc_UMAP_scRNA_cell_types.pdf')
DimPlot(object = pbmc_sample_umap, reduction = "umap",label = F,label.size = 4,repel=T,cols=c('B'='#1AB4B8','T'='#ED877F','mac'="#7DC05B",
'DC'='#9E9F20','NK'='#E06AA4'))+ NoLegend()
dev.off()

sample_cell<-Idents(pbmc_sample_umap)
cell_summary<-data.frame(strtrim(names(sample_cell),7),sample_cell)
colnames(cell_summary)<-c('sample_name','cell_type')
save(cell_summary,file='cell_summary.Rdata')

pdf('feature_gene_dot.pdf',width=10,height=10)
pbmc_sample_umap@active.ident<-factor(pbmc_sample_umap@active.ident,levels=rev(c('T','B','NK','mac','DC')))
DotPlot(pbmc_sample_umap,features=rev(feature_gene),cols=c('#BEBEBE','#FF0000'))
dev.off()

table(pbmc_sample_umap)
