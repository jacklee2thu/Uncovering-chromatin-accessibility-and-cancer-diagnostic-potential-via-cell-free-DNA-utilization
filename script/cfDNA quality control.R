samtools='/fshare2/lijie/download_software/samtools_1.15.1/bin/samtools'
cd /Sshare/lijie/cfDNA_LJ_protocol/datasets/quality_control
${samtools} flagstat WangLiXue_sort.bam > /Sshare/lijie/cfDNA_LJ_protocol/result/quality_control/WangLiXue_sort_flag.txt
${samtools} flagstat WangLiXue.bam > /Sshare/lijie/cfDNA_LJ_protocol/result/quality_control/WangLiXue_flag.txt
${samtools} coverage WangLiXue.bam > /Sshare/lijie/cfDNA_LJ_protocol/result/quality_control/WangLiXue_coverage.txt

options(stringsAsFactors=F)
setwd('E:/cfDNA_LJ_protocol/result/quality_control')
matrix_reads<-data.frame()
LJ_10_12_he_coverage<-read.table(file="LJ_10_12_he_coverage.txt",sep="\t",header=T)
LJ_10_12_he_sort_flag<-read.table(file="LJ_10_12_he_sort_flag.txt",sep="\t",header=F)
LJ_10_12_he_flag<-read.table(file="LJ_10_12_he_flag.txt",sep="\t",header=F)
LJ_10_12_he_sort_flagnum<-as.numeric(lapply(strsplit(LJ_10_12_he_sort_flag[,1],'+',fixed = T),function(x){return(x[[1]])}))
LJ_10_12_he_flagnum<-as.numeric(lapply(strsplit(LJ_10_12_he_flag[,1],'+',fixed = T),function(x){return(x[[1]])}))
sum_reads<-c(LJ_10_12_he_sort_flagnum[1],LJ_10_12_he_flagnum[1],LJ_10_12_he_flagnum[5])
matrix_reads<-rbind(matrix_reads,sum_reads)
colnames(matrix_reads)<-c('all reads','all deduplicated reads','mapped genome reads')
rownames(matrix_reads)<-c('GuJian','WangLiXue','LJ_10_12_he')
save(matrix_reads,GuJian_coverage,LJ_10_12_he_coverage,WangLiXue_coverage,file='matrix_final.Rdata')

options(stringsAsFactors=F)
setwd('E:/cfDNA_LJ_protocol/result/quality_control')
load('matrix_final.Rdata')
dfm<-data.frame()
for(i in 1:3){
temp_dfm<-data.frame(rownames(matrix_reads),matrix_reads[,i],rep(colnames(matrix_reads)[i],3))
dfm<-rbind(dfm,temp_dfm)
}
colnames(dfm)<-c('Description','reads','group')
dfm[,3]<-factor(dfm[,3],levels=c('all reads','all deduplicated reads','mapped genome reads'))

library(ggplot2)
library(ggpubr)
library(ggrepel)

ggbarplot(dfm, x = "Description", y = "reads", 
          add = c("mean_se"),sort.by.groups=T,
          color = "group",fill = "group", palette = "npg",
position = position_dodge(0.8))

LJ_10_12_he_coverage<-LJ_10_12_he_coverage[1:25,]
ggbarplot(LJ_10_12_he_coverage, x = "rname", y = "numreads", 
          add = c("mean_se"),sort.by.groups=T,
          color = "red2",fill = "red2", palette = "npg",
position = position_dodge(0.8))
ggbarplot(LJ_10_12_he_coverage, x = "rname", y = "meanmapq", 
          add = c("mean_se"),sort.by.groups=T,
          color = "red2",fill = "red2", palette = "npg",
position = position_dodge(0.8))
ggbarplot(LJ_10_12_he_coverage, x = "rname", y = "coverage", 
          add = c("mean_se"),sort.by.groups=T,
          color = "red2",fill = "red2", palette = "npg",
position = position_dodge(0.8))
ggbarplot(LJ_10_12_he_coverage, x = "rname", y = "meandepth", 
          add = c("mean_se"),sort.by.groups=T,
          color = "red2",fill = "red2", palette = "npg",
position = position_dodge(0.8))
