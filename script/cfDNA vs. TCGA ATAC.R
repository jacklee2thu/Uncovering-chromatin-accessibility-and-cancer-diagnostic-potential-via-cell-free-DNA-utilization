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

raw_insert<-readRDS('/Share2/home/lanxun3/jacklee/cfDNA/PANcancer_ATAC/TCGA-ATAC_PanCan_Log2Norm_Counts.rds')
all_bed<-read.table(file='/Share2/home/lanxun3/jacklee/cfDNA/PANcancer_ATAC/tss_pro_table_select.bed',header=F)
signature_liftover<-read.table(file='/Share2/home/lanxun3/jacklee/cfDNA/PANcancer_ATAC/tss_pro_table_liftover.bed',sep='\t',header=F)
colnames(all_bed)<-c(colnames(raw_insert[,1:7]),c('chromosome','start','end','strand','index','gene_name'))
promoter_bed<-all_bed[all_bed$annotation%in%'Promoter',]

# load('/Share2/home/lanxun3/jacklee/cfDNA/gbm/STAD_signature_genes.Rdata')
promoter_bed<-promoter_bed[promoter_bed$gene_name%in%signature_liftover$V6,]

promoter_bed<-rbind(promoter_bed,all_bed[all_bed$gene_name%in%setdiff(signature_liftover$V6,unique(promoter_bed$gene_name)),])
index<-unlist(lapply(strsplit(colnames(raw_insert)[8:803],'_'),function(x){x[[1]]}))
cancer_raw_insert<-raw_insert[match(promoter_bed$name,raw_insert$name),c(8:803)]

exp_gene_index<-match(exp_gene,promoter_bed$gene_name)
STAD_exp_gene_raw_insert<-cancer_raw_insert[exp_gene_index,grep('STAD',index)]
mean_STAD_exp_gene<-apply(STAD_exp_gene_raw_insert,1,mean)

unexp_gene_index<-match(unexp_gene,promoter_bed$gene_name)
STAD_unexp_gene_raw_insert<-cancer_raw_insert[unexp_gene_index,grep('STAD',index)]
mean_STAD_unexp_gene<-apply(STAD_unexp_gene_raw_insert,1,mean)
wilcox.test(mean_STAD_exp_gene,mean_STAD_unexp_gene)

df_fc<-data.frame(ATAC=c(mean_STAD_exp_gene,mean_STAD_unexp_gene),type=factor(c(rep('mean_STAD_exp_gene',length(mean_STAD_exp_gene)),
rep('mean_STAD_unexp_gene',length(mean_STAD_unexp_gene))),levels=c('mean_STAD_exp_gene','mean_STAD_unexp_gene')))

library(ggplot2)
library(ggpubr)
color = c("#EE3536","#3A429B")
pdf('TCGA_STAD_ATAC_tss.pdf')
my_comparisons <- list( c('mean_STAD_exp_gene','mean_STAD_unexp_gene') )
df_fc %>% ggplot(aes(x=type,y=ATAC,fill=type))+geom_violin()+
  geom_boxplot(width=0.2,position=position_dodge(0.9),outlier.colour = NA,fill="white")+
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
  stat_compare_means(label.y = 7)
dev.off()
