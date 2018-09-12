library(DESeq2)
library(limma)
library(AnnotationDbi)
library(org.Hs.eg.db)
#read the data and bind different conditions into the same table

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
