library(edgeR)
library(limma)
rna_seq<-read.table("GSE60424_GEOSubmit_FC1to11_normalized_counts.txt",header =TRUE)
column_names<-c("whole_blood","Neutrophil","Monocytes","Bcells","CD4T","CD8T","NK")
healthy_1<-cbind(rna_seq$lib221,rna_seq[,7:12])
healthy_2<-cbind(rna_seq[,74],rna_seq[,68:73])
healthy_3<-cbind(rna_seq[,81],rna_seq[,75:80])
healthy_4<-cbind(rna_seq[,95],rna_seq[,89:94])
sepsis_1<-cbind(rna_seq[,40],rna_seq[,34:39])
sepsis_2<-cbind(rna_seq[,46],rna_seq[,41:45])
sepsis_3<-cbind(rna_seq[,53],rna_seq[,47:52])
summerized<-cbind(healthy_1,healthy_2,healthy_3,healthy_4)
summerized<-cbind(summerized,sepsis_1,sepsis_2,sepsis_3)
colnames(sepsis_2)<-column_names

colnames(summerized)<- c("healthy1,whole_blood","healthy1,Neutrophil","healthy1,Monocytes","healthy1,Bcells","healthy1,CD4T",
                       "healthy1,CD8T","healthy1,NK","healthy2,whole_blood","healthy2,Neutrophil","healthy2,Monocytes","healthy2,Bcells","healthy2,CD4T",
                       "healthy2,CD8T","healthy2,NK","healthy3,whole_blood","healthy3,Neutrophil","healthy3,Monocytes","healthy3,Bcells","healthy3,CD4T",
                       "healthy3,CD8T","healthy3,NK","healthy4,whole_blood","healthy4,Neutrophil","healthy4,Monocytes","healthy4,Bcells","healthy4,CD4T",
                       "healthy4,CD8T","healthy4,NK","sepsis1,whole_blood","sepsis1,Neutrophil","sepsis1,Monocytes","sepsis1,Bcells","sepsis1,CD4T",
                       "sepsis1,CD8T","sepsis1,NK","sepsis2,whole_blood","sepsis2,Neutrophil","sepsis2,Monocytes","sepsis2,Bcells","sepsis2,CD4T",
                       "sepsis2,CD8T","sepsis3,whole_blood","sepsis3,Neutrophil","sepsis3,Monocytes","sepsis3,Bcells","sepsis3,CD4T",
                       "sepsis3,CD8T","sepsis3,NK")
row.names(summerized)<-rna_seq[,1]
#names(info) <- c("Donor","Component")
component<-factor(c("whole_blood","Neutrophil","Monocytes","Bcells","CD4T",
                    "CD8T","NK","whole_blood","Neutrophil","Monocytes","Bcells",
                    "CD4T","CD8T","NK","whole_blood","Neutrophil","Monocytes","Bcells",
                    "CD4T","CD8T","NK","whole_blood","Neutrophil","Monocytes","Bcells",
                    "CD4T","CD8T","NK"))
donor <-factor(c("healthy_1","healthy_1","healthy_1","healthy_1","healthy_1","healthy_1","healthy_1",
                 "healthy_2","healthy_2","healthy_2","healthy_2","healthy_2","healthy_2","healthy_2",
                 "healthy_3","healthy_3","healthy_3","healthy_3","healthy_3","healthy_3","healthy_3",
                 "healthy_4","healthy_4","healthy_4","healthy_4","healthy_4","healthy_4","healthy_4"))
design<-model.matrix(~0 + component + donor)



y<-DGEList(summerized,group = component)
row.names(y$counts)<-rna_seq[,1]
y <- calcNormFactors(y)
y <- estimateDisp(y, design, robust=TRUE)
fit <- glmQLFit(y, design, robust=TRUE)
#y$common.dispersion
qlf <- glmQLFTest(fit)
topTags(qlf,n=150)

# based on different individual
x<-DGEList(summerized,group = donor)
row.names(x$counts)<-rna_seq[,1]
x <- calcNormFactors(x)
x <- estimateDisp(x, design, robust=TRUE)
fit_x <- glmQLFit(x, design, robust=TRUE)
#y$common.dispersion
qlf <- glmQLFTest(fit)
topTags(qlf,n=150)

voom()
