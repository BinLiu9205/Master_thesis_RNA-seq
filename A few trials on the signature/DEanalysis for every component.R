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

whole_blood<-summerized[,grep("CD8T",colnames(summerized))]
group<-factor(c(rep("sepsis",3),rep("healthy",4)))

design<-model.matrix(~0 + group)
y<-DGEList(whole_blood,group = group)

 
 
 result<-voom(y,design,plot = T)
 
 fit <- lmFit(result, design)
 cont.matrix <- makeContrasts(Compare=groupsepsis-grouphealthy, levels=design)
 fit2 <- contrasts.fit(fit, cont.matrix)
 
 fit <- eBayes(fit2)
 # result<-topTable(fit, coef=ncol(design),number = 50045,sort.by = "p")
 
 result<-topTable(fit, number = 50045,sort.by = "p")
 result_sig<- result[result$adj.P.Val<0.05,]
 
 result_sig$genesymbol<-mapIds(org.Hs.eg.db, row.names(result_sig), column = "SYMBOL", "ENSEMBL")
 write.table(result,file="sepsis_vs_healthy_CD8T.txt",quote = FALSE) 
 write.table(as.data.frame(result_sig),file="sepsis_vs_healthy_CD8T_significant.txt",quote = FALSE)
 
 