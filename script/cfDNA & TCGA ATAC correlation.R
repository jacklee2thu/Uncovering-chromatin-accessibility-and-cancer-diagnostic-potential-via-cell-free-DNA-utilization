options(stringsAsFactors=F)
setwd('/Share2/home/lanxun3/jacklee/cfDNA/PANcancer_ATAC')
load('/Share2/home/lanxun3/jacklee/cfDNA/gbm/gastric_sample_tss_150.Rdata')
gastric_tss_150<-temp_tss_150
# load('/Share2/home/lanxun3/jacklee/cfDNA/gbm/gastric_sample_tss_1000.Rdata')
# gastric_tss_150<-temp_tss_1000
gastric_tss_150[is.na(gastric_tss_150)]<-0
load('/Share2/home/lanxun3/jacklee/pancancer_cfDNA/gbm/TGY_sample_tss_150.Rdata')
BRCA_tss_150<-temp_tss_150[,1:11]
# load('/Share2/home/lanxun3/jacklee/pancancer_cfDNA/gbm/TGY_sample_tss_1000.Rdata')
# BRCA_tss_150<-temp_tss_1000[,1:11]
BRCA_tss_150[is.na(BRCA_tss_150)]<-0
load('/Share2/home/lanxun3/jacklee/cfDNA/gbm/healthy_sample_tss_150.Rdata')
healthy_tss_150<-temp_tss_150
healthy_tss_150[is.na(healthy_tss_150)]<-0

load("/Share2/home/lanxun3/jacklee/cfDNA/tss/all_pro_tss/tss_pro_table.Rdata")
temp_tss<-tss_pro_table
temp_bed<-temp_tss[,c(1,2,3,5,8)]
temp_bed<-data.frame(temp_bed,index=1:70641)
write.table(temp_bed,file='tss_pro_table.bed',sep='\t',row.names=F,col.names=T,quote=F)
#http://genome.ucsc.edu/cgi-bin/hgLiftOver
signature_liftover<-read.table(file='tss_pro_table_liftover.bed',sep='\t',header=F)

write.table(raw_insert[,1:7],file='pancancer_bed.bed',sep='\t',row.names=F,col.names=F,quote=F)
bedtools intersect -a pancancer_bed.bed -b tss_pro_table_liftover.bed -wo >tss_pro_table_select.bed

#raw_insert<-readRDS('TCGA-ATAC_PanCan_Raw_Counts.rds')
raw_insert<-readRDS('TCGA-ATAC_PanCan_Log2Norm_Counts.rds')
all_bed<-read.table(file='tss_pro_table_select.bed',header=F)
colnames(all_bed)<-c(colnames(raw_insert[,1:7]),c('chromosome','start','end','strand','index','gene_name'))
promoter_bed<-all_bed[all_bed$annotation%in%'Promoter',]

# load('/Share2/home/lanxun3/jacklee/cfDNA/gbm/STAD_signature_genes.Rdata')
promoter_bed<-promoter_bed[promoter_bed$gene_name%in%signature_liftover$V6,]

promoter_bed<-rbind(promoter_bed,all_bed[all_bed$gene_name%in%setdiff(signature_liftover$V6,unique(promoter_bed$gene_name)),])
index<-unlist(lapply(strsplit(colnames(raw_insert)[8:803],'_'),function(x){x[[1]]}))
cancer_raw_insert<-raw_insert[match(promoter_bed$name,raw_insert$name),c(8:803)]

STAD_raw_insert<-cancer_raw_insert[,grep('STAD',index)]
BRCA_raw_insert<-cancer_raw_insert[,grep('BRCA',index)]
mean_STAD<-apply(STAD_raw_insert,1,mean)
mean_BRCA<-apply(BRCA_raw_insert,1,mean)
blood_STAD<-apply(gastric_tss_150[promoter_bed$index,],1,mean)
blood_BRCA<-apply(BRCA_tss_150[promoter_bed$index,],1,mean)
blood_healthy<-apply(healthy_tss_150[promoter_bed$index,],1,mean)
promoter_bed$mean_STAD<-mean_STAD
promoter_bed$blood_STAD<-blood_STAD
promoter_bed$mean_BRCA<-mean_BRCA
promoter_bed$blood_BRCA<-blood_BRCA
promoter_bed$blood_healthy<-blood_healthy
STAD_promoter_frame<-data.frame()
for(i in 1:length(unique(promoter_bed$gene_name))){
temp_promoter<-promoter_bed[promoter_bed$gene_name%in%unique(promoter_bed$gene_name)[i],]
temp_promoter<-temp_promoter[temp_promoter$mean_STAD==max(temp_promoter$mean_STAD),]
temp_promoter<-temp_promoter[temp_promoter$blood_STAD==min(temp_promoter$blood_STAD),]
STAD_promoter_frame<-rbind(STAD_promoter_frame,temp_promoter)
}
BRCA_promoter_frame<-data.frame()
for(i in 1:length(unique(promoter_bed$gene_name))){
temp_promoter<-promoter_bed[promoter_bed$gene_name%in%unique(promoter_bed$gene_name)[i],]
temp_promoter<-temp_promoter[temp_promoter$mean_BRCA==max(temp_promoter$mean_BRCA),]
temp_promoter<-temp_promoter[temp_promoter$blood_BRCA==min(temp_promoter$blood_BRCA),]
BRCA_promoter_frame<-rbind(BRCA_promoter_frame,temp_promoter)
}
save(STAD_promoter_frame,BRCA_promoter_frame,file='promoter_frame.Rdata')
cor.test(STAD_promoter_frame$blood_STAD,STAD_promoter_frame$mean_STAD)
cor.test(BRCA_promoter_frame$blood_BRCA,BRCA_promoter_frame$mean_BRCA)
cor.test(STAD_promoter_frame$blood_healthy,STAD_promoter_frame$mean_STAD)
cor.test(BRCA_promoter_frame$blood_healthy,BRCA_promoter_frame$mean_BRCA)



library(ggplot2)
dfm<-data.frame(new_blood_BRCA=c(STAD_promoter_frame$blood_STAD,BRCA_promoter_frame$blood_BRCA),
new_mean_BRCA=c(STAD_promoter_frame$mean_STAD,BRCA_promoter_frame$mean_BRCA),group=c(rep('STAD',length(STAD_promoter_frame$mean_STAD)),
rep('BRCA',length(BRCA_promoter_frame$mean_BRCA))))
pdf('blood_ATAC_allgene_cor.pdf',width=11,height=10)
p<-ggplot(dfm,aes(x=new_blood_BRCA,y=new_mean_BRCA,fill=group,color=group))+geom_smooth(method='lm')+geom_point(size=0.5,alpha=0.3)
p+scale_color_manual(values=c("#EE3536","#3A429B"))+theme_bw()+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())
dev.off()
