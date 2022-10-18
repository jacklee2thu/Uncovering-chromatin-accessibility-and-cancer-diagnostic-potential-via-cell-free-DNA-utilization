options(stringsAsFactors=F)
setwd('/Share2/home/lanxun3/jacklee/cfDNA/tss/all_pro_tss')
load("/Share2/home/lanxun3/jacklee/cfDNA/tss/all_pro_tss/tss_pro_table.Rdata")
load('/Share2/home/lanxun3/jacklee/cfDNA/tss/all_pro_tss/gastric_sample_tss_150.Rdata')
load('/Share2/home/lanxun3/jacklee/cfDNA/tss/all_pro_tss/gastric_sample_tss_1000.Rdata')

fre_150_1<-apply(temp_tss_150,1,function(x){length(which(x<1|is.na(x)==T))/dim(temp_tss_150)[2]})
names(fre_150_1)<-tss_pro_table$gene_name
fre_1000_1<-apply(temp_tss_1000,1,function(x){length(which(x<1|is.na(x)==T))/dim(temp_tss_1000)[2]})
names(fre_1000_1)<-tss_pro_table$gene_name
exp_gene<-intersect(names(which(fre_1000_1>=0.8)),names(which(fre_150_1>=0.8)))

setwd('/Share2/home/lanxun3/jacklee/cfDNA/tss/all_pro_tss')
library(pheatmap)
pdf('tss1000_relative_coverage.pdf')
new_tss_1000<-temp_tss_1000[order(fre_1000_1,decreasing = T),]
new_tss_1000[new_tss_1000>=2]<-2
pheatmap(new_tss_1000,cluster_cols = F,cluster_rows = F,
color = colorRampPalette(c("#D45252","white", "#0000B2"))(50),show_colnames = F)
dev.off()
pdf('tss150_relative_coverage.pdf')
new_tss_150<-temp_tss_150[order(fre_150_1,decreasing = T),]
new_tss_150[new_tss_150>=2]<-2
pheatmap(new_tss_150,cluster_cols = F,cluster_rows = F,
color = colorRampPalette(c("#D45252","white", "#0000B2"))(50),show_colnames = F)
dev.off()

index_sample_cancer_1000<-rep(0,length(fre_150_1))
index_sample_cancer_150<-rep(0,length(fre_150_1))
index_sample_cancer_150[match(intersect(which(fre_1000_1>=0.8),which(fre_150_1>=0.8)),order(fre_150_1,decreasing = T))]<-1
index_sample_cancer_1000[match(intersect(which(fre_1000_1>=0.8),which(fre_150_1>=0.8)),order(fre_1000_1,decreasing = T))]<-1

gene_set<-unique(exp_gene)
library(org.Hs.eg.db)
k=keys(org.Hs.eg.db,keytype = "ENSEMBL")
list=select(org.Hs.eg.db,keys=k,columns = c("ENTREZID","SYMBOL"), keytype="ENSEMBL")
most_123456_genename=na.omit(list[match(gene_set,list[,"SYMBOL"]),][,2])
library(clusterProfiler)
KEGG_gene<-enrichKEGG(most_123456_genename,organism = "human", pvalueCutoff = 0.05)
GO_gene<-enrichGO(most_123456_genename,OrgDb=org.Hs.eg.db,ont = "BP", pvalueCutoff = 0.05)
write.table(KEGG_gene@result,file='KEGG_tss.txt',sep='\t',row.names = F,col.names = T,quote = F)
write.table(GO_gene@result,file='GO_tss.txt',sep='\t',row.names = F,col.names = T,quote = F)


View(KEGG_gene@result)
View(GO_gene@result)

options(stringsAsFactors=F)
setwd('E:/gastric cfDNA/results/TSS/tss_coverage/case_figure')
KEGG_tss<-read.table(file='KEGG_tss.txt',sep='\t',header=T,quote='')
GO_tss<-read.table(file='GO_tss.txt',sep='\t',header=T,quote='')
observed_KEGG<-c('Gastric cancer','Epstein-Barr virus infection','Shigellosis',
'Epithelial cell signaling in Helicobacter pylori infection','mTOR signaling pathway','Cell cycle','EGFR tyrosine kinase inhibitor resistance',
'ErbB signaling pathway','AMPK signaling pathway','HIF-1 signaling pathway','MAPK signaling pathway',
'Wnt signaling pathway','TNF signaling pathway','VEGF signaling pathway','TGF-beta signaling pathway',
'Endocytosis','T cell receptor signaling pathway')

all_path<-GO_tss
path_index<-match(observed_GO,all_path[,2])
new_kegg<-all_path[path_index,c(2,5,9)]
all_kegg<-new_kegg

final_frame<-all_kegg
final_frame<-na.omit(final_frame)
#final_frame$Description<-factor(final_frame$Description,levels=rev(observed_GO))
final_frame$Description<-factor(final_frame$Description,levels=final_frame[order(final_frame$Count),]$Description)
library(ggplot2)
p = ggplot(final_frame,aes(Count,Description))
p = ggplot(final_frame,aes(Count,Description))
p=p + geom_point()  
p=p + geom_point(aes(size=Count))
pbubble = p+ geom_point(aes(size=Count,color=pvalue))
pr = pbubble+scale_color_gradient(low="red",high = "yellow")
pr = pr+labs(color=expression(pvalue),size="Count",  
                           x="GeneRatio",y="Pathway name",title="Pathway enrichment")
pr + theme_bw()

pdf('KEGG_exp.pdf')
dotplot(KEGG_gene,showCategory=30)
dev.off()
pdf('GO_exp.pdf',width=15,height=10)
dotplot(GO_gene,showCategory=30)
dev.off()

####
fre_150<-apply(temp_tss_150,1,function(x){length(which(x>=1.5))/dim(temp_tss_150)[2]})
names(fre_150)<-tss_pro_table$gene_name
fre_1000<-apply(temp_tss_1000,1,function(x){length(which(x>=1.5))/dim(temp_tss_1000)[2]})
names(fre_1000)<-tss_pro_table$gene_name
unexp_gene<-intersect(names(which(fre_1000>=0.1)),names(which(fre_150>=0.1)))#954
unexp_gene<-setdiff(unexp_gene,intersect(exp_gene,unexp_gene))

index_sample_cancer_150[match(intersect(which(fre_1000>=0.1),which(fre_150>=0.1)),order(fre_150_1,decreasing = T))]<-2
index_sample_cancer_1000[match(intersect(which(fre_1000>=0.1),which(fre_150>=0.1)),order(fre_1000_1,decreasing = T))]<-2

index_sample<-data.frame(index_sample_cancer_150,index_sample_cancer_1000)
pdf('index_sample.pdf')
pheatmap(index_sample,cluster_cols = F,cluster_rows = F,show_colnames = F,show_rownames = F,
color = colorRampPalette(c("#BEBEBE","#D45252",'#FFA500'))(50))
dev.off()


TCGA_count<-read.table(file="/Share2/home/lanxun3/jacklee/cfDNA/tss/all_pro_tss/HTSeq - Counts.merge.txt",sep="\t",header = T)
rownames(TCGA_count)<-TCGA_count[,1]
rownames(TCGA_count)<-strtrim(rownames(TCGA_count),15)
TCGA_count<-TCGA_count[,-1]
TCGA_patient<-TCGA_count[,colnames(TCGA_count)[which(substr(colnames(TCGA_count),14,15)=="01")]]

gene_set<-unique(exp_gene)
ungene_set<-unique(unexp_gene)
library(org.Hs.eg.db)
k=keys(org.Hs.eg.db,keytype = "ENSEMBL")
list=select(org.Hs.eg.db,keys=k,columns = c("ENTREZID","SYMBOL"), keytype="ENSEMBL")
gene_set_ensemble=na.omit(list[match(gene_set,list[,"SYMBOL"]),][,1])
ungene_set_ensemble=na.omit(list[match(ungene_set,list[,"SYMBOL"]),][,1])

gene_exp_TCGA<-TCGA_patient[gene_set_ensemble,]
ungene_exp_TCGA<-TCGA_patient[ungene_set_ensemble,]
gene_exp_mean<-apply(gene_exp_TCGA,1,mean)
ungene_exp_mean<-apply(ungene_exp_TCGA,1,mean)
wilcox.test(gene_exp_mean,ungene_exp_mean,alternative='great')
                                                                                                                                               
df<-data.frame(c(gene_exp_mean,ungene_exp_mean),c(rep("pre_exp",length(gene_exp_mean)),
rep("pre_unexp",length(ungene_exp_mean))))
colnames(df)<-c("num","group")

library(ggpubr)

library(digest)
setwd('/Share2/home/lanxun3/jacklee/cfDNA/tss/all_pro_tss')
pdf('tss_TCGA_exp.pdf')
p<-ggboxplot(df, "group", "num",color = "group", size = 1,width=0.7,outlier.colour=NA, palette =c("#D45252",'#FFA500',"#6EA6D9"),add = "jitter", shape = "group")
p+scale_y_log10()
dev.off()
