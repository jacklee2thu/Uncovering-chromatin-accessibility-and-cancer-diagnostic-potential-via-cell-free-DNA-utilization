options(stringsAsFactors=F)
setwd("E:/gastric cfDNA/results/TSS/TCGA_tss")
load('E:/gastric cfDNA/GTEx/blood_GTEX_file_count.Rdata')
load("tss_pro_table.Rdata")
TCGA_count<-read.table(file="E:/gastric cfDNA/TCGA_count/count/HTSeq - Counts.merge.txt",sep="\t",header = T)
#TCGA_count<-read.table(file="E:/breast cancer/TCGA/TCGA_RNA_seq_Count.txt",sep="\t",header = T)
rownames(TCGA_count)<-TCGA_count[,1]
rownames(TCGA_count)<-strtrim(rownames(TCGA_count),15)
TCGA_count<-TCGA_count[,-1]
pcg<-read.table("E:/科研/reference_genome/hg_19_pro.txt",sep = "\t",header = T)
TCGA_count<-TCGA_count[intersect(rownames(TCGA_count),pcg$Gene.stable.ID),]
Gtex_blood<-blood_GTEX_file_count[intersect(rownames(blood_GTEX_file_count),pcg$Gene.stable.ID),]
TCGA_patient<-TCGA_count[,colnames(TCGA_count)[which(substr(colnames(TCGA_count),14,15)=="01")]]
TCGA_normal<-TCGA_count[,colnames(TCGA_count)[which(substr(colnames(TCGA_count),14,15)=="11")]]
patient_count<-apply(TCGA_patient,1,mean)
patient_count<-patient_count[intersect(names(patient_count),tss_pro_table$ensemble_id)]

up_genes<-names(sort(patient_count,decreasing = T)[1:500])
down_genes<-names(sort(patient_count,decreasing = T)[18089:18588])
save(up_genes,down_genes,file='TCGA_exp.Rdata')


##Gtex_blood<-Gtex_blood[intersect(rownames(Gtex_blood),rownames(TCGA_patient)),]
##TCGA_patient<-TCGA_patient[intersect(rownames(Gtex_blood),rownames(TCGA_patient)),]
library(limma)
group_list=c(rep(1,dim(TCGA_patient)[2]),rep(0,dim(TCGA_normal)[2]))
exprSet<-cbind(TCGA_patient,TCGA_normal)
design <- model.matrix(~factor(group_list))
colnames(design)=levels(factor(group_list))
rownames(design)=colnames(exprSet)
v <- voom(exprSet,design,normalize="quantile")
fit <- lmFit(v,design)
fit2 <- eBayes(fit)
tempOutput = topTable(fit2, coef=2, n=Inf)
DEG_voom = na.omit(tempOutput)
save(DEG_voom,file="TCGA_diff.Rdata")
#save(DEG_voom,file="TCGA_BRCA_diff.Rdata")
new_signature<-c('RAE1','GPR85','CALU','AK9','TRAF2','PAPD4','ABCA10','LRRIQ1','SPECC1','MEIS2')
new_signature<-c('CACNA1B','ASAP2','KREMEN1','FGFR3','CDH4','SRCIN1','SNRPN','CTD-2021H9.3','ANKS1B','WNT4')

#######火山图
temp_frame<-DEG_voom
fc=2
temp_frame$threshold<-0
temp_frame$threshold[temp_frame$logFC > fc & temp_frame$adj.P.Val <= 0.01] = "up"
temp_frame$threshold[temp_frame$logFC < -fc & temp_frame$adj.P.Val <= 0.01] = "down"
temp_frame$id<-rownames(DEG_voom)
temp_frame$id<-pcg[match(rownames(temp_frame),pcg$Gene.stable.ID),]$Gene.name


library(ggplot2)
library(ggrepel)
gg<-ggplot(data = temp_frame, aes(x = logFC, y = -log10(adj.P.Val),color=threshold)) +
  geom_point(size=3) +
  geom_vline(xintercept=c(-fc,fc),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(0.01),lty=4,col="black",lwd=0.8) +
  labs(x="log2_FC",y="-log10_adj.P.Val",title="TCGA_RNA_deg") +
  theme(plot.title = element_text(hjust = 0.5), legend.position="right", legend.title = element_blank()) +  
    geom_text_repel(
    data = temp_frame[temp_frame$p_value_fdr>=2.193226,],
    aes(label = id),
    size = 5,
    box.padding = unit(0.5, "lines"),
    point.padding = unit(0.8, "lines"), segment.color = "black", show.legend = FALSE )
gg+theme_bw()+theme(axis.line = element_line(size=1, colour = "black"),panel.border = element_blank()
,panel.grid.major=element_blank(),panel.grid.minor=element_blank())


###TSS
options(stringsAsFactors=F)
setwd("/Share2/home/lanxun3/jacklee/cfDNA/tss/tcga_tss")
load("tss_pro_table.Rdata")
#load("TCGA_diff.Rdata")
# load('TCGA_BRCA_diff.Rdata')
#load("TCGA_tumor_Gtex_blood.Rdata")
load("TCGA_exp.Rdata")
# up_genes<-rownames(DEG_voom[DEG_voom$adj.P.Val<=0.01&DEG_voom$logFC>=2,])
# down_genes<-rownames(DEG_voom[DEG_voom$adj.P.Val<=0.01&DEG_voom$logFC<= -2,])

TCGA_tss_up<-c()#####2856 5359
for(i in 1:length(up_genes)){
temp_index<-which(tss_pro_table$ensemble_id==up_genes[i])
TCGA_tss_up<-c(TCGA_tss_up,temp_index)
}
TCGA_tss_down<-c()#####596 3610
for(i in 1:length(down_genes)){
temp_index<-which(tss_pro_table$ensemble_id==down_genes[i])
TCGA_tss_down<-c(TCGA_tss_down,temp_index)
}
TCGA_tss_index<-c(TCGA_tss_up,TCGA_tss_down)


setwd("/Share2/home/lanxun3/jacklee/cfDNA/tss/all_pro_tss")
load('gastric_sample_tss_150.Rdata')
load('gastric_sample_tss_1000.Rdata')
temp_tss_150<-temp_tss_150[TCGA_tss_index,]
temp_tss_1000<-temp_tss_1000[TCGA_tss_index,]

save(temp_tss_150,temp_tss_1000,file="/Share2/home/lanxun3/jacklee/cfDNA/tss/tcga_tss/tcga_tss.Rdata")
#save(temp_tss_150,temp_tss_1000,file="/Share2/home/lanxun3/jacklee/cfDNA/tss/tcga_tss/tcga_diff_tss.Rdata")
####boxplot
setwd("E:/gastric cfDNA/results/TSS/TCGA_tss")
load("tcga_tss.Rdata")
# load("tcga_BRCA_diff_tss.Rdata")
#load("tcga_diff_tss.Rdata")
temp_tss_150[is.na(temp_tss_150)]<-0.1
temp_tss_1000[is.na(temp_tss_1000)]<-0.1
high100_1000<-temp_tss_1000[1:2856,]
low100_1000<-temp_tss_1000[2856:3452,]
#high100_1000<-temp_tss_1000[1:1431,]
#low100_1000<-temp_tss_1000[1432:3347,]
a<-apply(high100_1000,1,mean)
b<-apply(low100_1000,1,mean)
wilcox.test(a,b,alternative = 'less')

high100_150<-temp_tss_150[1:2856,]
low100_150<-temp_tss_150[2856:3452,]
#high100_150<-temp_tss_150[1:1431,]
#low100_150<-temp_tss_150[1432:3347,]

a<-apply(high100_150,1,mean)
b<-apply(low100_150,1,mean)
wilcox.test(a,b,alternative = 'less')

df<-data.frame(c(a,b),c(rep("high500",length(a)),
rep("low500",length(b))))
colnames(df)<-c("num","group")

library(ggpubr)

library(digest)
p<-ggboxplot(df, "group", "num",color = "group", size = 1,width=0.7,outlier.colour=NA, palette =c( "#D45252",'#FFA500',"#6EA6D9"),add = "jitter", shape = "group")
p+ labs(title = 'coverage')+stat_compare_means()+ylim(0,2)
