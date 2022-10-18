argv <- commandArgs(TRUE)
sample_name_prefix <- as.character(argv[1])
options(stringsAsFactors=F)
##
setwd("/fshare2/lijie/PAN_CANCER_cfDNA/datasets")
patient_sample<-list.dirs(full.names = FALSE,recursive = F)
#for(k in 1:length(patient_sample)){
k<-grep(sample_name_prefix,patient_sample)
setwd(paste0("/fshare2/lijie/PAN_CANCER_cfDNA/datasets/",patient_sample[k]))
if(is.na(match(paste0(patient_sample[k],".bam.bai"),dir()))){
system(paste0("/Share2/home/lanxun/Tools/samtools-1.3.1/bin/bin/samtools index ",patient_sample[k],".bam"))
}

tss_bed<-read.table(file="/Share2/home/lanxun3/jacklee/cfDNA/tss/immune/tss_pro_bed_up3000_1000.bed",sep="\t",header=F)
tss_bed_up3000_1000<-read.table(file="/Share2/home/lanxun3/jacklee/cfDNA/tss/immune/tss_pro_bed_up3000_1000.bed",sep="\t",header=F)
tss_bed_down1000_3000<-read.table(file="/Share2/home/lanxun3/jacklee/cfDNA/tss/immune/tss_pro_bed_down1000_3000.bed",sep="\t",header=F)
tss_bed_up1000_down1000<-read.table(file="/Share2/home/lanxun3/jacklee/cfDNA/tss/immune/tss_pro_bed_up1000_down1000.bed",sep="\t",header=F)
tss_bed_up150_down50<-read.table(file="/Share2/home/lanxun3/jacklee/cfDNA/tss/immune/tss_pro_bed_up150_down50.bed",sep="\t",header=F)

patient_tss1000<-list()
patient_tss150<-list()
for(i in 1:dim(tss_bed)[1]){
##
temp_file<-system(paste0("/Share2/home/lanxun/Tools/samtools-1.3.1/bin/bin/samtools",
" depth -q 30"," -r ",tss_bed_up3000_1000[i,1],":",tss_bed_up3000_1000[i,2],"-",tss_bed_up3000_1000[i,3]," ",patient_sample[k],".bam"),intern=TRUE, wait=TRUE)

if(length(temp_file)==0){
temp_up3000_1000<-NA
}else{
temp_up3000_1000<-data.frame()
for(j in 1:length(temp_file)){
temp_frame<-unlist(strsplit(temp_file[j],'\t'))
temp_up3000_1000<-rbind(temp_up3000_1000,temp_frame)
}
colnames(temp_up3000_1000)<-c('chr','site','coverage')
temp_up3000_1000[,3]<-as.numeric(temp_up3000_1000[,3])
}

temp_file<-system(paste0("/Share2/home/lanxun/Tools/samtools-1.3.1/bin/bin/samtools",
" depth -q 30"," -r ",tss_bed_down1000_3000[i,1],":",tss_bed_down1000_3000[i,2],"-",tss_bed_down1000_3000[i,3]," ",patient_sample[k],".bam"),intern=TRUE, wait=TRUE)
if(length(temp_file)==0){
temp_down1000_3000<-NA
}else{
temp_down1000_3000<-data.frame()
for(j in 1:length(temp_file)){
temp_frame<-unlist(strsplit(temp_file[j],'\t'))
temp_down1000_3000<-rbind(temp_down1000_3000,temp_frame)
}
colnames(temp_down1000_3000)<-c('chr','site','coverage')
temp_down1000_3000[,3]<-as.numeric(temp_down1000_3000[,3])
}

if(unique(is.na(temp_down1000_3000))&unique(is.na(temp_up3000_1000))){
normal_tss_value<-NA
}
if(unique(is.na(temp_up3000_1000))&unique(!is.na(temp_down1000_3000))){
normal_tss_value<-mean(temp_down1000_3000[,3])
}
if(unique(is.na(temp_down1000_3000))&unique(!is.na(temp_up3000_1000))){
normal_tss_value<-mean(temp_up3000_1000[,3])
}

if(unique(!is.na(temp_down1000_3000))&unique(!is.na(temp_up3000_1000))){
normal_tss<-rbind(temp_up3000_1000,temp_down1000_3000)
normal_tss_value<-mean(normal_tss[,3])
}


##
temp_file<-system(paste0("/Share2/home/lanxun/Tools/samtools-1.3.1/bin/bin/samtools",
" depth -q 30"," -r ",tss_bed_up1000_down1000[i,1],":",tss_bed_up1000_down1000[i,2],"-",tss_bed_up1000_down1000[i,3]," ",patient_sample[k],".bam"),intern=TRUE, wait=TRUE)
if(length(temp_file)==0){
temp_up1000_down1000<-NA
}else{
temp_up1000_down1000<-data.frame()
for(j in 1:length(temp_file)){
temp_frame<-unlist(strsplit(temp_file[j],'\t'))
temp_up1000_down1000<-rbind(temp_up1000_down1000,temp_frame)
}
temp_up1000_down1000[,3]<-as.numeric(temp_up1000_down1000[,3])
temp_up1000_down1000[,4]<-temp_up1000_down1000[,3]/normal_tss_value
colnames(temp_up1000_down1000)<-c('chr','site','coverage','mean_coverage')
}

##
temp_file<-system(paste0("/Share2/home/lanxun/Tools/samtools-1.3.1/bin/bin/samtools",
" depth -q 30"," -r ",tss_bed_up150_down50[i,1],":",tss_bed_up150_down50[i,2],"-",tss_bed_up150_down50[i,3]," ",patient_sample[k],".bam"),intern=TRUE, wait=TRUE)
if(length(temp_file)==0){
temp_up150_down50<-NA
}else{
temp_up150_down50<-data.frame()
for(j in 1:length(temp_file)){
temp_frame<-unlist(strsplit(temp_file[j],'\t'))
temp_up150_down50<-rbind(temp_up150_down50,temp_frame)
}
temp_up150_down50[,3]<-as.numeric(temp_up150_down50[,3])
temp_up150_down50[,4]<-temp_up150_down50[,3]/normal_tss_value
colnames(temp_up150_down50)<-c('chr','site','coverage','mean_coverage')
}

patient_tss1000[[i]]<-temp_up1000_down1000
patient_tss150[[i]]<-temp_up150_down50
print(i)
}
save(patient_tss1000,patient_tss150,file=paste0(patient_sample[k],"_alltss.Rdata"))

options(stringsAsFactors=F)
setwd("/Share2/home/lanxun3/jacklee/cfDNA/gastric_cancer")
last_sample<-dir()
load("/Share2/home/lanxun3/jacklee/cfDNA/tss/all_pro_tss/tss_pro_table.Rdata")
all_sample_tss<-list()
for(j in 1:length(last_sample)){
load(paste0(last_sample[j],'/',last_sample[j],"_alltss.Rdata"))
tss_1000<-c()
tss_150<-c()

for(i in 1:length(patient_tss1000)){
if(!is.na(patient_tss1000[[i]])){
tss_1000[i]<-mean(patient_tss1000[[i]][,4])
}else{tss_1000[i]<-NA}
if(!is.na(patient_tss150[[i]])){
tss_150[i]<-mean(patient_tss150[[i]][,4])
}else{tss_150[i]<-NA}
}
tss_relative_coverage<-data.frame(tss_150,tss_1000)
colnames(tss_relative_coverage)<-c("tss_150","tss_1000")
all_sample_tss[[j]]<-tss_relative_coverage
}


names(all_sample_tss)<-last_sample
setwd('/Share2/home/lanxun3/jacklee/cfDNA/tss/all_pro_tss')
save(all_sample_tss,file='all_sample_tss.Rdata')
#save(all_sample_tss,file='healthy_sample_tss.Rdata')


###########
temp_tss_150<-all_sample_tss[[1]][,1]
temp_tss_1000<-all_sample_tss[[1]][,2]
for(i in 2:length(all_sample_tss)){
temp_tss_150<-cbind(temp_tss_150,all_sample_tss[[i]][,1])
temp_tss_1000<-cbind(temp_tss_1000,all_sample_tss[[i]][,2])
}
colnames(temp_tss_150)<-last_sample
colnames(temp_tss_1000)<-last_sample
save(temp_tss_150,file='gastric_sample_tss_150.Rdata')
save(temp_tss_1000,file='gastric_sample_tss_1000.Rdata')
#save(temp_tss_150,file='healthy_sample_tss_150.Rdata')
#save(temp_tss_1000,file='healthy_sample_tss_1000.Rdata')
