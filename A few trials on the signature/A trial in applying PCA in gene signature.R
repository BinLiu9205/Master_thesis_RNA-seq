rna_seq<-read.table("GSE60424_GEOSubmit_FC1to11_normalized_counts.txt",header =TRUE)
column_names<-c("whole_blood","Neutrophil","Monocytes","Bcells","CD4T","CD8T","NK")
healthy_1<-cbind(rna_seq$lib221,rna_seq[,7:12])
healthy_2<-cbind(rna_seq[,74],rna_seq[,68:73])
healthy_3<-cbind(rna_seq[,81],rna_seq[,75:80])
healthy_4<-cbind(rna_seq[,95],rna_seq[,89:94])
sepsis_1<-cbind(rna_seq[,40],rna_seq[,34:39])
sepsis_2<-cbind(rna_seq[,46],rna_seq[,41:45])
sepsis_3<-cbind(rna_seq[,53],rna_seq[,47:52])

summerized<-cbind(sepsis_1,sepsis_2,sepsis_3)
summerized<-cbind(summerized,healthy_1,healthy_2,healthy_3,healthy_4)

#colnames(sepsis_1)<-column_names 

colnames(summerized)<- c("sepsis1,whole_blood","sepsis1,Neutrophil","sepsis1,Monocytes","sepsis1,Bcells","sepsis1,CD4T",
                         "sepsis1,CD8T","sepsis1,NK","sepsis2,whole_blood","sepsis2,Neutrophil","sepsis2,Monocytes","sepsis2,Bcells","sepsis2,CD4T",
                         "sepsis2,CD8T","sepsis3,whole_blood","sepsis3,Neutrophil","sepsis3,Monocytes","sepsis3,Bcells","sepsis3,CD4T",
                         "sepsis3,CD8T","sepsis3,NK","healthy1,whole_blood","healthy1,Neutrophil","healthy1,Monocytes","healthy1,Bcells","healthy1,CD4T",
                         "healthy1,CD8T","healthy1,NK","healthy2,whole_blood","healthy2,Neutrophil","healthy2,Monocytes","healthy2,Bcells","healthy2,CD4T",
                         "healthy2,CD8T","healthy2,NK","healthy3,whole_blood","healthy3,Neutrophil","healthy3,Monocytes","healthy3,Bcells","healthy3,CD4T",
                         "healthy3,CD8T","healthy3,NK","healthy4,whole_blood","healthy4,Neutrophil","healthy4,Monocytes","healthy4,Bcells","healthy4,CD4T",
                         "healthy4,CD8T","healthy4,NK")
row.names(summerized)<-rna_seq[,1]
result <- as.data.frame(t(log2(1+summerized)))
result$disease <- c("sepsis","sepsis","sepsis","sepsis","sepsis",
                    "sepsis","sepsis","sepsis","sepsis","sepsis","sepsis","sepsis",
                    "sepsis","sepsis","sepsis","sepsis","sepsis","sepsis",
                    "sepsis","sepsis","healthy","healthy","healthy","healthy","healthy",
                    "healthy","healthy","healthy","healthy","healthy","healthy","healthy",
                    "healthy","healthy","healthy","healthy","healthy","healthy","healthy",
                    "healthy","healthy","healthy","healthy","healthy","healthy","healthy",
                    "healthy","healthy")
result$cell_type <-  c("whole_blood","Neutrophil","Monocytes","Bcells","CD4T",
                       "CD8T","NK","whole_blood","Neutrophil","Monocytes","Bcells","CD4T",
                       "CD8T","whole_blood","Neutrophil","Monocytes","Bcells","CD4T",
                       "CD8T","NK","whole_blood","Neutrophil","Monocytes","Bcells","CD4T",
                       "CD8T","NK","whole_blood","Neutrophil","Monocytes","Bcells","CD4T",
                       "CD8T","NK","whole_blood","Neutrophil","Monocytes","Bcells","CD4T",
                       "CD8T","NK","whole_blood","Neutrophil","Monocytes","Bcells","CD4T",
                       "CD8T","NK")
PCA_all<- prcomp(result[,1:50045])
library(ggplot2)
ggplot(data=NULL, aes(PCA_all$x[,1],PCA_all$x[,2])) + geom_point(aes(color=result$cell_type,shape=result$disease))+
  xlab("First Dimension") + ylab("Second Dimension") + labs(shape="Health_Status",color="Cell_type")
PoV <- PCA_all$sdev^2/sum(PCA_all$sdev^2)*100
barplot(PoV, xlab= "Dimensions", ylab="Proportion of explained variance (%)") 
PoV[1]+PoV[2]
######The healthy candidates only, but with partial components
result_healthy_PC <- result[result$disease=="healthy"&result$cell_type!="whole_blood",]
PCA_healthy <- prcomp(result_healthy_PC[,1:50045])
ggplot(data=NULL, aes(PCA_healthy$x[,1],PCA_healthy$x[,2])) + geom_point(aes(color=result_healthy_PC$cell_type))+
  xlab("First Dimension") + ylab("Second Dimension") + labs(color="Cell_type")
PoV <- PCA_healthy$sdev^2/sum(PCA_healthy$sdev^2)*100
barplot(PoV, xlab= "Dimensions", ylab="Proportion of explained variance (%)") 
PoV[1]+PoV[2]
target_gene <- result[result$disease=="healthy"&result$cell_type=="whole_blood",]
score_x <- NULL
score_y <-NULL
for (i in 1:4){
  score_x[i] <- (as.numeric(target_gene[i,1:50045]) - PCA_healthy$center) %*% PCA_healthy$rotation[,1]
  score_y[i] <- (as.numeric(target_gene[i,1:50045]) - PCA_healthy$center) %*% PCA_healthy$rotation[,2] 
}
mean(score_x)
mean(score_y)


result_sepsis_PC <- result[result$disease=="sepsis"&result$cell_type!="whole_blood",]
PCA_sepsis <- prcomp(result_sepsis_PC[,1:50045])
ggplot(data=NULL, aes(PCA_sepsis$x[,1],PCA_sepsis$x[,2])) + geom_point(aes(color=result_sepsis_PC$cell_type))+
  xlab("First Dimension") + ylab("Second Dimension") + labs(color="Cell_type")
PoV <- PCA_sepsis$sdev^2/sum(PCA_sepsis$sdev^2)*100
barplot(PoV, xlab= "Dimensions", ylab="Proportion of explained variance (%)") 
PoV[1]+PoV[2]
target_s_gene <- result[result$disease=="sepsis"&result$cell_type=="whole_blood",]
score_s_x <- NULL
score_s_y <-NULL
for (i in 1:3){
  score_s_x[i] <- (as.numeric(target_s_gene[i,1:50045]) - PCA_sepsis$center) %*% PCA_sepsis$rotation[,1]
  score_s_y[i] <- (as.numeric(target_s_gene[i,1:50045]) - PCA_sepsis$center) %*% PCA_sepsis$rotation[,2] 
}
mean(score_s_x)
mean(score_s_y)

healthy_res <- as.data.frame(cbind(score_x,score_y))
sepsis_res <- as.data.frame(cbind(score_s_x,score_s_y))
colnames(healthy_res)<- c("PC1","PC2")
colnames(sepsis_res)<- c("PC1","PC2")
all_res <- as.data.frame(rbind(healthy_res,sepsis_res))
all_res$health <- c(rep("healthy",4),rep("sepsis",3))
ggplot(data = all_res, mapping = aes(PC1, PC2))+geom_point(aes(color=all_res$health))+
  labs(color="healthy status")
############################################################

#Using different components in blood to analyze the result
result_healthy_PC <- result[result$disease=="healthy"&(result$cell_type=="Monocytes"|result$cell_type=="Bcells"|result$cell_type=="Neutrophil"),]
PCA_healthy <- prcomp(result_healthy_PC[,1:50045])
ggplot(data=NULL, aes(PCA_healthy$x[,1],PCA_healthy$x[,2])) + geom_point(aes(color=result_healthy_PC$cell_type))+
  xlab("First Dimension") + ylab("Second Dimension") + labs(color="Cell_type")
PoV <- PCA_healthy$sdev^2/sum(PCA_healthy$sdev^2)*100
barplot(PoV, xlab= "Dimensions", ylab="Proportion of explained variance (%)") 
PoV[1]+PoV[2]
target_gene <- result[result$disease=="healthy"&result$cell_type=="whole_blood",]
score_x <- NULL
score_y <-NULL
for (i in 1:4){
  score_x[i] <- (as.numeric(target_gene[i,1:50045]) - PCA_healthy$center) %*% PCA_healthy$rotation[,1]
  score_y[i] <- (as.numeric(target_gene[i,1:50045]) - PCA_healthy$center) %*% PCA_healthy$rotation[,2] 
}
mean(score_x)
mean(score_y)

#result$cell_type=="Neutrophil"|
result_sepsis_PC <- result[result$disease=="sepsis"&(result$cell_type=="Monocytes"|result$cell_type=="Bcells"|result$cell_type=="Neutrophil"),]
PCA_sepsis <- prcomp(result_sepsis_PC[,1:50045])
ggplot(data=NULL, aes(PCA_sepsis$x[,1],PCA_sepsis$x[,2])) + geom_point(aes(color=result_sepsis_PC$cell_type))+
  xlab("First Dimension") + ylab("Second Dimension") + labs(color="Cell_type")
PoV <- PCA_sepsis$sdev^2/sum(PCA_sepsis$sdev^2)*100
barplot(PoV, xlab= "Dimensions", ylab="Proportion of explained variance (%)") 
PoV[1]+PoV[2]
target_s_gene <- result[result$disease=="sepsis"&result$cell_type=="whole_blood",]
score_s_x <- NULL
score_s_y <-NULL
#prediction process 
for (i in 1:3){
  score_s_x[i] <- (as.numeric(target_s_gene[i,1:50045]) - PCA_sepsis$center) %*% PCA_sepsis$rotation[,1]
  score_s_y[i] <- (as.numeric(target_s_gene[i,1:50045]) - PCA_sepsis$center) %*% PCA_sepsis$rotation[,2] 
}
mean(score_s_x)
mean(score_s_y)

healthy_res <- as.data.frame(cbind(score_x,score_y))
sepsis_res <- as.data.frame(cbind(score_s_x,score_s_y))
colnames(healthy_res)<- c("PC1","PC2")
colnames(sepsis_res)<- c("PC1","PC2")
all_res <- as.data.frame(rbind(healthy_res,sepsis_res))
all_res$health <- c(rep("healthy",4),rep("sepsis",3))
ggplot(data = all_res, mapping = aes(PC1, PC2))+geom_point(aes(color=all_res$health))+
  labs(color="healthy status")


#####To work on which component contribute the most to the classification###
W_blood <- result[result$cell_type=="whole_blood",]
s_W_blood <- W_blood[W_blood$disease=="sepsis",]
h_W_blood <- W_blood[W_blood$disease=="healthy",]
s_score <-as.data.frame(predict(PCA_all,s_W_blood)) 
h_score <- as.data.frame(predict(PCA_all,h_W_blood))
all_s <-as.data.frame(rbind(s_score[,1:2],h_score[,1:2])) 
all_s$type <- c(rep("sepsis",3),rep("healthy",4))
ggplot(all_s,aes(x=PC1,y=PC2)) +geom_point(aes(color=all_s$type))+
  labs(color="healthy status")
s_Neutr <- result[result$cell_type=="Neutrophil"&result$disease=="sepsis",] 
h_Neutr <- result[result$cell_type=="Neutrophil"&result$disease=="healthy",] 
s_neu_score <-as.data.frame(predict(PCA_all,s_Neutr)) 
h_neu_score <- as.data.frame(predict(PCA_all,h_Neutr))
all_Neutr <- as.data.frame((rbind(s_neu_score[,1:2],h_neu_score[,1:2])) )
all_Neutr$type <- c(rep("sepsis",3),rep("healthy",4))
s_Bcells <- result[result$cell_type=="Bcells"&result$disease=="sepsis",] 
h_Bcells <- result[result$cell_type=="Bcells"&result$disease=="healthy",] 
s_Bcells_score <-as.data.frame(predict(PCA_all,s_Bcells)) 
h_Bcells_score <- as.data.frame(predict(PCA_all,h_Bcells))
all_Bcells <- as.data.frame((rbind(s_Bcells_score[,1:2],h_Bcells_score[,1:2])) )
all_Bcells$type <- c(rep("sepsis",3),rep("healthy",4))
s_Monocytes <- result[result$cell_type=="Monocytes"&result$disease=="sepsis",] 
h_Monocytes <- result[result$cell_type=="Monocytes"&result$disease=="healthy",] 
s_Mono_score <-as.data.frame(predict(PCA_all,s_Monocytes)) 
h_Mono_score <- as.data.frame(predict(PCA_all,h_Monocytes))
all_Mono <- as.data.frame((rbind(s_Mono_score[,1:2],h_Mono_score[,1:2])) )
all_Mono$type <- c(rep("sepsis",3),rep("healthy",4))
over_all <- rbind(all_s,all_Neutr,all_Bcells,all_Mono)
over_all $ celltype <- c(rep("whole_blood",7),rep("neutrophil",7),rep("Bcells",7),rep("Monocytes",7))
ggplot(over_all,aes(x=PC1,y=PC2)) +geom_point(aes(color=over_all$type,shape=over_all$celltype))+
  labs(color="healthy status",shape="cell_type")
