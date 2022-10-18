options(stringsAsFactors=F)
setwd('/Share2/home/lanxun3/jacklee/cfDNA/gbm')
patient_info<-read.table(file='patient_info.txt',sep='\t',header=T,quote='',na.strings='')
load('gastric_sample_tss_150.Rdata')
gastric_tss_150<-temp_tss_150
gastric_sample<-colnames(gastric_tss_150)

signet_ring_sample<-patient_info[patient_info$histology%in%c('Signet-ring cell carcinoma','signet-ring cell carcinoma'),]$research
Adenocarcinoma_sample<-patient_info[patient_info$histology%in%c('Adenocarcinoma'),]$research
signet_ring_tss_150<-gastric_tss_150[,intersect(gastric_sample,c(signet_ring_sample,'XieXiangRong'))]
rownames(signet_ring_tss_150)<-paste0("tss_150_", 1:70641)
signet_ring_tss_150[is.na(signet_ring_tss_150)]<-0
signet_ring_tss_150<-as.data.frame(signet_ring_tss_150)
Adenocarcinoma_tss_150<-gastric_tss_150[,intersect(gastric_sample,Adenocarcinoma_sample)]
rownames(Adenocarcinoma_tss_150)<-paste0("tss_150_", 1:70641)
Adenocarcinoma_tss_150[is.na(Adenocarcinoma_tss_150)]<-0
Adenocarcinoma_tss_150<-as.data.frame(Adenocarcinoma_tss_150)


library(limma)
group_list=c(rep(1,dim(signet_ring_tss_150)[2]),rep(0,dim(Adenocarcinoma_tss_150)[2]))
exprSet<-cbind(signet_ring_tss_150,Adenocarcinoma_tss_150)
design <- model.matrix(~factor(group_list))
colnames(design)=levels(factor(group_list))
rownames(design)=colnames(exprSet)
v <- voom(exprSet,design,normalize="quantile")
fit <- lmFit(v,design)
fit2 <- eBayes(fit)
tempOutput = topTable(fit2, coef=2, n=Inf)
DEG_voom = na.omit(tempOutput)
save(DEG_voom,file='TCGA/signet_ring_vs_Adenocarcinoma_DEG.Rdata')
DEG_voom$signet_mean<-apply(signet_ring_tss_150[rownames(DEG_voom),],1,mean)
DEG_voom$Adenocarcinoma_mean<-apply(Adenocarcinoma_tss_150[rownames(DEG_voom),],1,mean)
######火山图
signet_ring_up_signature<-rownames(DEG_voom[DEG_voom$logFC<=0&DEG_voom$P.Value<=0.05&DEG_voom$signet_mean<1&DEG_voom$Adenocarcinoma_mean>1,])
Adenocarcinoma_up_signature<-rownames(DEG_voom[DEG_voom$logFC>=0&DEG_voom$P.Value<=0.05&DEG_voom$Adenocarcinoma_mean<1&DEG_voom$signet_mean>1,])
DEG_voom$threshold<-rep(0,dim(DEG_voom)[1])
DEG_voom[signet_ring_up_signature,]$threshold<-1
DEG_voom[Adenocarcinoma_up_signature,]$threshold<-2
DEG_voom$threshold<-as.factor(DEG_voom$threshold)
key_DEG_voom<-DEG_voom[c(signet_ring_up_signature[1:10],Adenocarcinoma_up_signature[1:10]),]

index<-unlist(lapply(strsplit(rownames(key_DEG_voom),'_'),function(x){x[[3]]}))
load("/Share2/home/lanxun3/jacklee/cfDNA/tss/all_pro_tss/tss_pro_table.Rdata")
signature_gene<-tss_pro_table[as.numeric(index),]$gene_name
df<-data.frame(key_DEG_voom$logFC,-log10(key_DEG_voom$P.Value),signature_gene)
colnames(df)<-c('log2FC','-log10(P.Value)','gene_name')
pdf('TCGA/volcano_signet_ring_vs_Adenocarcinoma_raw.pdf')
library(ggplot2)
library(ggrepel)
gg<-ggplot(DEG_voom, aes(x = -logFC, y = -log10(P.Value),color=threshold)) +
  geom_point(size=1) +
  geom_hline(yintercept = -log10(0.05),lty=4,col="black",lwd=0.8)+
  scale_color_manual(breaks = c(0,1,2),values=c("#BEBEBE","#D45252","#6EA6D9"))
gg+theme_bw()+theme(axis.line = element_line(size=1, colour = "black"),panel.border = element_blank()
,panel.grid.major=element_blank(),panel.grid.minor=element_blank(),legend.position="none")
dev.off()
pdf('TCGA/volcano_signet_ring_vs_Adenocarcinoma_dot.pdf')
library(ggplot2)
gg<-ggplot(df,aes(x=log2FC,y= df[,2]))+geom_point(size=1)+annotate("text",
x=df$log2FC,y= df[,2],label=df$gene_name)
gg<-gg+ geom_hline(aes(yintercept=0), colour="black")+geom_vline(aes(xintercept=0), colour="black")
gg+theme_bw()+theme(axis.line = element_line(size=1, colour = "black"),panel.border = element_blank()
,panel.grid.major=element_blank(),panel.grid.minor=element_blank(),legend.position="none")
dev.off()
##########功能富集
library(org.Hs.eg.db)
k=keys(org.Hs.eg.db,keytype = "ENSEMBL")
list=select(org.Hs.eg.db,keys=k,columns = c("ENTREZID","SYMBOL"), keytype="ENSEMBL")
####signet_ring_up_signature_function
index<-unlist(lapply(strsplit(signet_ring_up_signature,'_'),function(x){x[[3]]}))
load("/Share2/home/lanxun3/jacklee/cfDNA/tss/all_pro_tss/tss_pro_table.Rdata")
up_gastric_genes<-tss_pro_table[as.numeric(index),]$gene_name
pro_entrezid=list[match(unique(up_gastric_genes),list[,"SYMBOL"]),][,2]

library(clusterProfiler)
KEGG_gene<-enrichKEGG(pro_entrezid,organism = "hsa", pvalueCutoff = 0.05)
KEGG_gene<-setReadable(KEGG_gene,OrgDb='org.Hs.eg.db',keyType='ENTREZID')
GO_gene<-enrichGO(pro_entrezid,'org.Hs.eg.db',ont = "BP", pvalueCutoff = 0.05)
GO_gene<-setReadable(GO_gene,OrgDb='org.Hs.eg.db',keyType='ENTREZID')
setwd('/Share2/home/lanxun3/jacklee/cfDNA/gbm/TCGA')
write.table(KEGG_gene@result,file='signet_ring_up_signature_KEGG_gene.txt',sep='\t',row.names = F,col.names = T,quote = F)
write.table(GO_gene@result,file='signet_ring_up_signature_GO_gene.txt',sep='\t',row.names = F,col.names = T,quote = F) 
save(KEGG_gene,GO_gene,file='signet_ring_up_signature_function.Rdata')
######Adenocarcinoma_up_signature_function
index<-unlist(lapply(strsplit(Adenocarcinoma_up_signature,'_'),function(x){x[[3]]}))
load("/Share2/home/lanxun3/jacklee/cfDNA/tss/all_pro_tss/tss_pro_table.Rdata")
up_gastric_genes<-tss_pro_table[as.numeric(index),]$gene_name
pro_entrezid=list[match(unique(up_gastric_genes),list[,"SYMBOL"]),][,2]

library(clusterProfiler)
KEGG_gene<-enrichKEGG(pro_entrezid,organism = "hsa", pvalueCutoff = 0.05)
KEGG_gene<-setReadable(KEGG_gene,OrgDb='org.Hs.eg.db',keyType='ENTREZID')
GO_gene<-enrichGO(pro_entrezid,'org.Hs.eg.db',ont = "BP", pvalueCutoff = 0.05)
GO_gene<-setReadable(GO_gene,OrgDb='org.Hs.eg.db',keyType='ENTREZID')
setwd('/Share2/home/lanxun3/jacklee/cfDNA/gbm/TCGA')
write.table(KEGG_gene@result,file='Adenocarcinoma_up_signature_KEGG_gene.txt',sep='\t',row.names = F,col.names = T,quote = F)
write.table(GO_gene@result,file='Adenocarcinoma_up_signature_GO_gene.txt',sep='\t',row.names = F,col.names = T,quote = F) 
save(KEGG_gene,GO_gene,file='Adenocarcinoma_up_signature_function.Rdata')
#####气泡图
setwd('/Share2/home/lanxun3/jacklee/cfDNA/gbm/TCGA')
load('signet_ring_up_signature_function.Rdata')
load('Adenocarcinoma_up_signature_function.Rdata')
observed_path<-c('morphogenesis of a branching epithelium','branching morphogenesis of an epithelial tube',
'morphogenesis of a branching structure','cell-matrix adhesion','cell-substrate adhesion',
'positive regulation of cytokine production','regulation of leukocyte activation')

i=2
all_path<-rbind(GO_gene@result,KEGG_gene@result)

path_index<-match(observed_path,all_path[,2])
new_kegg<-data.frame(all_path[path_index,c(2,5,8,9)],group=rep(i,length(observed_path)))
#all_kegg<-new_kegg
all_kegg<-rbind(all_kegg,new_kegg)
save(all_kegg,file="signet_ring_all_kegg.Rdata")
final_frame<-all_kegg
final_frame<-na.omit(final_frame)
final_frame$pvalue[final_frame$pvalue> 0.1]<-0.1
final_frame$Description<-factor(final_frame$Description,levels=observed_path)
final_frame$group<-factor(final_frame$group,levels=c(2,1))


pdf('signet_ring_vs_Adenocarcinoma_pathway.pdf',width=12,height=6)
library(ggplot2)
p = ggplot(final_frame,aes(Description,group))
p=p + geom_point()  
# 修稿点的大小
p=p + geom_point(aes(size=Count))
# 展示三维数据
pbubble = p+ geom_point(aes(size=Count,color=pvalue))
# 设置渐变色
pr = pbubble+scale_color_gradient(low="red",high = "blue")
# 绘制p气泡图
pr = pr+labs(color=expression(pvalue),size="Count",  
                           x="Pathway name",y="group",title="Pathway enrichment")
pr + theme_bw()
dev.off()

##ESR1,ESR2,CDH1,ERBB2
index<-which(tss_pro_table$gene_name=='ESR1')
i=9
df<-data.frame(c(as.numeric(signet_ring_tss_150[index[i],]),as.numeric(Adenocarcinoma_tss_150[index[i],])),
c(rep("signet_ring",length(signet_ring_tss_150[index[i],])),
rep("Adenocarcinoma",length(Adenocarcinoma_tss_150[index[i],]))))
colnames(df)<-c("num","group")
library(ggpubr)
library(digest)
pdf('ESR1_signet_ring_vs_Adenocarcinoma.pdf',width=10,height=10)
p<-ggboxplot(df, "group", "num",color = "group", palette =c("#D45252",'#FFA500'),add = "jitter")
p+stat_compare_means()
dev.off()
wilcox.test(as.numeric(signet_ring_tss_150[index[i],]),as.numeric(Adenocarcinoma_tss_150[index[i],]),alternative='less')
