options(stringsAsFactors=F)
setwd('/Share2/home/lanxun3/jacklee/cfDNA/experiment/ATAC/SGC7901_ATAC_2/peaks')
peak_bed<-read.table(file="SGC7901_ATAC_2_peak_peaks.narrowPeak",sep="\t",header=F)
colnames(peak_bed)<-c('chr','start','end','peak','-10*log10(q)','strand','FC','-log10(p)','-log10(q)','relative_summit')
peak_bed$up_peak<-peak_bed$start-1000
peak_bed$down_peak<-peak_bed$end+1000

peak_bed_up1000_start<-peak_bed[,c('chr','up_peak','start')]
peak_bed_start_end<-peak_bed[,c('chr','start','end')]
peak_bed_end_down1000<-peak_bed[,c('chr','end','down_peak')]

write.table(peak_bed_up1000_start,file="peak_bed_up1000_start.bed",sep="\t",quote = F,col.names =F,row.names=F)
write.table(peak_bed_start_end,file="peak_bed_start_end.bed",sep="\t",quote = F,col.names =F,row.names=F)
write.table(peak_bed_end_down1000,file="peak_bed_end_down1000.bed",sep="\t",quote = F,col.names =F,row.names=F)

##
options(stringsAsFactors=F)
setwd("/Share2/home/lanxun3/jacklee/cfDNA/experiment/cfDNA")
patient_sample<-list.dirs(full.names = FALSE,recursive = F)
#for(k in 1:length(patient_sample)){
k=1
setwd(paste0("/Share2/home/lanxun3/jacklee/cfDNA/experiment/cfDNA/",patient_sample[k]))
system(paste0("/Share2/home/lanxun/Tools/samtools-1.3.1/bin/bin/samtools index ",patient_sample[k],".bam"))
options(stringsAsFactors=F)
peak_bed<-read.table(file="/Share2/home/lanxun3/jacklee/cfDNA/experiment/ATAC/GSE_1_ATAC/peaks/peak_bed_up1000_start.bed",
sep="\t",header=F)

patient_peak<-list()
for(i in 1:dim(peak_bed)[1]){
##
peak_bed<-read.table(file="/Share2/home/lanxun3/jacklee/cfDNA/experiment/ATAC/GSE_1_ATAC/peaks/peak_bed_up1000_start.bed",
sep="\t",header=F)
system(paste0("/Share2/home/lanxun/Tools/samtools-1.3.1/bin/bin/samtools",
" depth -q 30"," -r ",peak_bed[i,1],":",peak_bed[i,2],"-",peak_bed[i,3]," ",patient_sample[k],".bam"," > ",patient_sample[k],"_temppeak.txt"),intern=TRUE, wait=TRUE)
if(file.info(paste0(patient_sample[k],"_temppeak.txt"))$size == 0){
temp_up1000_start<-NA
}else{
temp_up1000_start<-read.table(file=paste0(patient_sample[k],"_temppeak.txt"),sep="\t",head=F)
}
system(paste0("rm /Share2/home/lanxun3/jacklee/cfDNA/experiment/cfDNA/",
patient_sample[k],"/",patient_sample[k],"_temppeak.txt"))

peak_bed<-read.table(file="/Share2/home/lanxun3/jacklee/cfDNA/experiment/ATAC/GSE_1_ATAC/peaks/peak_bed_end_down1000.bed",sep="\t",header=F)
system(paste0("/Share2/home/lanxun/Tools/samtools-1.3.1/bin/bin/samtools",
" depth -q 30"," -r ",peak_bed[i,1],":",peak_bed[i,2],"-",peak_bed[i,3]," ",patient_sample[k],".bam"," > ",patient_sample[k],"_temppeak.txt"),intern=TRUE, wait=TRUE)
if(file.info(paste0(patient_sample[k],"_temppeak.txt"))$size == 0){
temp_end_down1000<-NA
}else{
temp_end_down1000<-read.table(file=paste0(patient_sample[k],"_temppeak.txt"),sep="\t",head=F)
}
system(paste0("rm /Share2/home/lanxun3/jacklee/cfDNA/experiment/cfDNA/",
patient_sample[k],"/",patient_sample[k],"_temppeak.txt"))
if(unique(is.na(temp_up1000_start))&unique(is.na(temp_end_down1000))){
normal_peak_value<-NA
}
if(unique(is.na(temp_up1000_start))&unique(!is.na(temp_end_down1000))){
normal_peak_value<-mean(temp_end_down1000[,3])
}
if(unique(!is.na(temp_up1000_start))&unique(is.na(temp_end_down1000))){
normal_peak_value<-mean(temp_up1000_start[,3])
}

if(unique(!is.na(temp_up1000_start))&unique(!is.na(temp_end_down1000))){
normal_peak<-rbind(temp_up1000_start,temp_end_down1000)
normal_peak_value<-mean(normal_peak[,3])
}


##start-end
peak_bed<-read.table(file="/Share2/home/lanxun3/jacklee/cfDNA/experiment/ATAC/GSE_1_ATAC/peaks/peak_bed_start_end.bed",sep="\t",header=F)
system(paste0("/Share2/home/lanxun/Tools/samtools-1.3.1/bin/bin/samtools",
" depth -q 30"," -r ",peak_bed[i,1],":",peak_bed[i,2],"-",peak_bed[i,3]," ",patient_sample[k],".bam"," > ",patient_sample[k],"_temppeak.txt"),intern=TRUE, wait=TRUE)
if(file.info(paste0(patient_sample[k],"_temppeak.txt"))$size == 0){
temp_start_end<-NA
}else{
temp_start_end<-read.table(file=paste0(patient_sample[k],"_temppeak.txt"),sep="\t",head=F)
temp_start_end[,4]<-temp_start_end[,3]/normal_peak_value
}
system(paste0("rm /Share2/home/lanxun3/jacklee/cfDNA/experiment/cfDNA/",
patient_sample[k],"/",patient_sample[k],"_temppeak.txt"))

patient_peak[[i]]<-temp_start_end
}
save(patient_peak,file=paste0(patient_sample[k],"_allpeak.Rdata"))

options(stringsAsFactors=F)
setwd('/Share2/home/lanxun3/jacklee/cfDNA/experiment/ATAC/SGC7901_ATAC_2/peaks')
peak_bed<-read.table(file="SGC7901_ATAC_2_peak_peaks.narrowPeak",sep="\t",header=F)
colnames(peak_bed)<-c('chr','start','end','peak','-10*log10(q)','strand','FC','-log10(p)','-log10(q)','relative_summit')
setwd('/Share2/home/lanxun3/jacklee/cfDNA/experiment/cfDNA/SGC7901_cfDNA_5h')
load('SGC7901_cfDNA_5h_allpeak.Rdata')
peak_bed$peak_location<-peak_bed$start+peak_bed$relative_summit
mean_raw<-c()
mean_relative<-c()
new_peak_100<-list()
for(i in 1:length(patient_peak)){
temp_peak<-patient_peak[[i]]
if(is.na(temp_peak)){
mean_raw[i]<-NA
mean_relative[i]<-NA
new_peak_100[[i]]<-NA
}else{
index1<-which(temp_peak[,2]>=peak_bed[i,]$peak_location-100)
index2<-which(temp_peak[,2]<=peak_bed[i,]$peak_location+100)
if(length(index1)!=0&length(index2)!=0){
new_peak<-temp_peak[min(index1):max(index2),]
new_peak_100[[i]]<-new_peak
new_peak_100[[i]]$V2<-new_peak_100[[i]]$V2-peak_bed[i,]$peak_location
mean_raw[i]<-mean(new_peak[,3])
mean_relative[i]<-mean(new_peak[,4])
}else{
mean_raw[i]<-NA
mean_relative[i]<-NA
new_peak_100[[i]]<-NA
}
}
}
save(new_peak_100,file='new_peak_100bp.Rdata',version =2)

###
options(stringsAsFactors=F)
setwd('/Share2/home/lanxun3/jacklee/cfDNA/experiment/cfDNA/MKN28_cfDNA_5h')
peak_bed<-read.table(file="/Share2/home/lanxun3/jacklee/cfDNA/experiment/ATAC/SGC7901_ATAC_2/peaks/SGC7901_ATAC_2_peak_peaks.narrowPeak",sep="\t",header=F)
colnames(peak_bed)<-c('chr','start','end','peak','-10*log10(q)','strand','FC','-log10(p)','-log10(q)','relative_summit')
peak_bed$peak_location<-peak_bed$start+peak_bed$relative_summit
write.table(peak_bed,file='peak_bed.bed',sep='\t',row.names=F,col.names=F,quote=F)
#GSE_1_ATAC_file<-read.table(file='/Share2/home/lanxun3/jacklee/cfDNA/experiment/ATAC/GSE_1_ATAC/peaks/GSE_1_ATAC_peak_treat_pileup.bdg',sep='\t',header=F)
format(Sys.time(), "%D:%H:%M:%S")
temp_coverage<-system(paste0("/Share2/home/lanxun/Tools/samtools-1.3.1/bin/bin/samtools",
" depth -q 30"," -b ","peak_bed.bed"," ",
"/Share2/home/lanxun3/jacklee/cfDNA/experiment/ATAC/SGC7901_ATAC_2/SGC7901_ATAC_2.last.bam"),intern=TRUE, wait=TRUE)
format(Sys.time(), "%D:%H:%M:%S")
GSE_1_ATAC_file<-data.frame()
for(j in 1:length(temp_coverage)){
temp_frame<-unlist(strsplit(temp_coverage[j],'\t'))
GSE_1_ATAC_file<-rbind(GSE_1_ATAC_file,temp_frame)
}
GSE_1_ATAC_file[,2]<-as.numeric(GSE_1_ATAC_file[,2])
GSE_1_ATAC_file[,3]<-as.numeric(GSE_1_ATAC_file[,3])
colnames(GSE_1_ATAC_file)<-c('chr','site','coverage')
format(Sys.time(), "%D:%H:%M:%S")
new_ATAC_100<-list()
for(i in 1:dim(peak_bed)[1]){
target_location<-peak_bed[i,]
temp_file<-GSE_1_ATAC_file[GSE_1_ATAC_file[,1]%in%target_location[,1],]
index1<-which(temp_file[,2]>=peak_bed[i,]$peak_location-100)
index2<-which(temp_file[,2]<=peak_bed[i,]$peak_location+100)
if(length(index1)!=0&length(index2)!=0){
new_frame<-temp_file[min(index1):max(index2),]
new_ATAC_100[[i]]<-new_frame
}else{
new_ATAC_100[[i]]<-NA
}
new_ATAC_100[[i]][,2]<-new_ATAC_100[[i]][,2]-peak_bed[i,]$peak_location
}
save(new_ATAC_100,file='new_ATAC_100bp.Rdata')
format(Sys.time(), "%D:%H:%M:%S")
new_ATAC_100_dataframe<-data.frame()
for(i in 1:length(new_ATAC_100)){
if(!is.na(new_ATAC_100[[i]])){
new_ATAC_100_dataframe<-rbind(new_ATAC_100_dataframe,new_ATAC_100[[i]])
}else{
new_ATAC_100_dataframe<-new_ATAC_100_dataframe
}
}
save(new_ATAC_100_dataframe,file='new_ATAC_100_dataframe.Rdata')
df<-data.frame(new_ATAC_100_dataframe$V2,new_ATAC_100_dataframe$V4)
colnames(df)<-c('location','fitted')
mean_relative_cfDNA<-c()
for(i in 1:201){
mean_relative_cfDNA[i]<-sum(df[df$location%in%(i-101),]$fitted)/length(new_ATAC_100)
}
new_df<-data.frame(c(-100:100),mean_relative_cfDNA)
colnames(new_df)<-c('location','mean_relative')
library(ggplot2)

pdf('MKN28_100bp_ATAC_loess.pdf')
p<-ggplot(new_df,aes(x=location,y=mean_relative))+geom_smooth(method = 'loess')+geom_vline(aes(xintercept=0),colour="#6EA6D9", linetype="dashed",size=1)
p+theme_bw()+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())
dev.off()
pdf('MKN28_100bp_ATAC_loess_1.pdf')
p<-ggplot(new_df,aes(x=location,y=mean_relative))+geom_line()+geom_vline(aes(xintercept=0),colour="#6EA6D9", linetype="dashed",size=1)
p+theme_bw()+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())
dev.off()
