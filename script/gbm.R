options(stringsAsFactors=F)
setwd('/Share2/home/lanxun3/jacklee/cfDNA/gbm')
patient_info<-read.table(file='patient_info.txt',sep='\t',header=T,quote='',na.strings='')
load('gastric_sample_tss_150.Rdata')
gastric_tss_150<-temp_tss_150
load('healthy_sample_tss_150.Rdata')
healthy_tss_150<-temp_tss_150
set.seed(1)
healthy_tss_150<-healthy_tss_150[,sample(colnames(healthy_tss_150),dim(healthy_tss_150)[2])]
load('NGY_sample_tss_150.Rdata')
NGY_tss_150<-temp_tss_150
load('TT_sample_tss_150.Rdata')
TT_tss_150<-temp_tss_150
load('GGY_sample_tss_150.Rdata')
GGY_tss_150<-temp_tss_150
features.sl <- rbind(t(gastric_tss_150),t(healthy_tss_150),t(NGY_tss_150),t(TT_tss_150),t(GGY_tss_150))
colnames(features.sl)<-paste0("tss_150_", 1:70641)
features.sl[is.na(features.sl)]<-0
features.sl<-as.data.frame(features.sl)
gastric_sample<-colnames(gastric_tss_150)
xiehe_sample<-colnames(cbind(healthy_tss_150,NGY_tss_150,TT_tss_150,GGY_tss_150))
#######stage
stage_sample<-patient_info[patient_info$Stage%in%c('Ia','Ib'),]$research
stage_sample<-patient_info[patient_info$Stage%in%c('IIa','IIb'),]$research
stage_sample<-patient_info[patient_info$Stage%in%c('IIIa','IIIb','IIIc','IV'),]$research
gastric_sample<-intersect(gastric_sample,stage_sample)
features.sl<-features.sl[c(gastric_sample,xiehe_sample),]
######gbm

library(tidyverse)
library(caret)
library(pROC)
library(gbm)
xiehe_index<-match(xiehe_sample,rownames(features.sl))
gastric_index<-match(gastric_sample,rownames(features.sl))
type_name<-rep(0,dim(features.sl)[1])
type_name[xiehe_index]<-'Healthy'
type_name[gastric_index]<-'Cancer'
type_name<-factor(type_name,levels=c('Cancer','Healthy'))
i<-12345693
set.seed(i)
inTrain = createDataPartition(type_name, p = 3/4, list = FALSE)
trainx = features.sl[inTrain,]
testx = features.sl[-inTrain,]
trainy = type_name[inTrain]
testy = type_name[-inTrain]

####

subsets = c(20,50,100)####
ctrl= rfeControl(functions = rfFuncs,method = 'cv',verbose = FALSE,returnResamp = 'final')
Profile = rfe(trainx, trainy, sizes = subsets, rfeControl = ctrl)
if(length(Profile$optVariables)>=101){
trainx=trainx[,Profile$optVariables[1:100]]
testx=testx[,Profile$optVariables[1:100]]
}else{
trainx=trainx[,Profile$optVariables]
testx=testx[,Profile$optVariables]
}
features<-trainx

features$type<-trainy

ctrl <- trainControl(method = "repeatedcv",
                     number = 10,
                     repeats = 10,
                     verboseIter = FALSE,
                     savePredictions=TRUE,
                     classProbs=TRUE,
#                      preProcOptions=list(thres = 0.90),
                     summaryFunction = twoClassSummary)
					 
gbmGrid = expand.grid(interaction.depth = 3,n.trees = 150,shrinkage = 0.1,n.minobsinnode = 10)
#gbmFit1 = train(features,method = 'gbm',trControl = ctrl,tuneGrid = gbmGrid,verbose = FALSE)
	
					 
					 
set.seed(1234)
model_gbm <- caret::train(type ~ .,
                               data = features,
                               method = 'gbm',
							   tuneGrid=gbmGrid,
							   preProcess = c("corr", "nzv"),
                         trControl = ctrl)



#### Save
save(model_gbm,file='model_gbm_tss150.Rdata')
pred.tbl <- model_gbm$pred %>%
group_by(rowIndex) %>% dplyr::summarize(obs=obs[1], Cancer=mean(Cancer))
pred.tbl$sample <- rownames(features)
combine_info<-cbind(rownames(features),features)
colnames(combine_info)[1]<-'sample'
pred.tbl <- inner_join(pred.tbl, combine_info)

## 90% specificity
cutoff <- (pred.tbl %>% filter(type=="Healthy") %>%
           arrange(desc(Cancer)))$Cancer[15]
## 95% specificity cutoff to be used in tissue prediction.
cutoff95 <- (pred.tbl %>% filter(type=="Healthy") %>%
             arrange(desc(Cancer)))$Cancer[9]

pred.tbl <- pred.tbl %>%
    mutate(detected90 = ifelse(Cancer > cutoff, "Cancer", "Healthy"),
       detected95 = ifelse(Cancer > cutoff95, "Cancer", "Healthy"))
	   
confusion_matrix<-confusionMatrix(factor(pred.tbl$detected90,levels=c('Cancer','Healthy')), pred.tbl$obs)
confusion_matrix
confusion_matrix<-confusionMatrix(factor(pred.tbl$detected95,levels=c('Cancer','Healthy')), pred.tbl$obs)
confusion_matrix
stage_sample<-patient_info[patient_info$Stage%in%c('Ia','Ib'),]$research
stage_sample<-patient_info[patient_info$Stage%in%c('IIa','IIb'),]$research
stage_sample<-patient_info[patient_info$Stage%in%c('IIIa','IIIb','IIIc','IV'),]$research
stage_pred.tbl<-pred.tbl[pred.tbl$sample%in%stage_sample,]
confusion_matrix<-confusionMatrix(factor(stage_pred.tbl$detected90,levels=c('Cancer','Healthy')), stage_pred.tbl$obs)
confusion_matrix
confusion_matrix<-confusionMatrix(factor(stage_pred.tbl$detected95,levels=c('Cancer','Healthy')), stage_pred.tbl$obs)
confusion_matrix
###
library(ROCR)
pred.tbl$lable=ifelse(pred.tbl$obs=='Cancer',yes=1,0)
pred1 = prediction(pred.tbl$Cancer,pred.tbl$lable)
perf1 = performance(pred1, measure='tpr', x.measure='fpr')
auc_ROCR <-performance(pred1,measure ="auc")
auc_ROCR@y.values[[1]]
save(perf1,file='ROC_cnv.Rdata')
library(ggplot2)
df <- data.frame(Curve=as.factor(rep(c(1), each=length(perf1@x.values[[1]]))), 
                 FalsePositive=c(perf1@x.values[[1]]),
                 TruePositive=c(perf1@y.values[[1]]))
pdf('tss150_STAD_train_ROC.pdf')
ggplot(df, aes(x=FalsePositive, y=TruePositive, color=Curve)) +ggtitle(paste0("ROC Curve  AUC=", auc_ROCR@y.values[[1]]))+ geom_line()+theme_bw()+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())
dev.off()
######

predict(model_gbm, newdata = testx)[1:5]
models<-list(model_gbm)
predValues = extractPrediction(models,testX = testx, testY = testy)
testValues = subset(predValues, dataType == 'Test')
probValues = extractProb(models,testX = testx, testY = testy)
testProbs = subset(probValues, dataType == 'Test')
pred.tbl = subset(testValues, model == 'gbm')
pred.tbl$sample<-rownames(testx)
#new_pred.tbl<-pred.tbl[c(grep('GGY',pred.tbl$sample),grep('TT',pred.tbl$sample),grep('Cancer',pred.tbl$obs)),]
confusion_matrix<-confusionMatrix(pred.tbl$pred, pred.tbl$obs)
confusion_matrix
save(confusion_matrix,file='confusion_matrix_tss150.Rdata')
library(ROCR)
pred.tbl$lable=ifelse(pred.tbl$obs=='Cancer',yes=1,0)
pred1 = prediction(testProbs$Cancer,pred.tbl$lable)
perf1 = performance(pred1, measure='tpr', x.measure='fpr')
auc_ROCR <-performance(pred1,measure ="auc")
auc_ROCR@y.values[[1]]
cbind(pred.tbl,rownames(testx))[,c(1,2,7)]
save(perf1,file='ROC_tss150.Rdata')

library(ggplot2)
df <- data.frame(Curve=as.factor(rep(c(1), each=length(perf1@x.values[[1]]))), 
                 FalsePositive=c(perf1@x.values[[1]]),
                 TruePositive=c(perf1@y.values[[1]]))
pdf('tss150_ROC.pdf')
ggplot(df, aes(x=FalsePositive, y=TruePositive, color=Curve)) + geom_line()+theme_bw()+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),legend.position='none')
dev.off()
##
temp_frame<-cbind(rownames(testx),pred.tbl)
early_frame<-rbind(temp_frame[temp_frame[,1]%in%intersect(temp_frame[,1],c(stage_sample1,stage_sample2)),],temp_frame[26:77,])
confusion_matrix<-confusionMatrix(early_frame[,3], early_frame[,2])
confusion_matrix
