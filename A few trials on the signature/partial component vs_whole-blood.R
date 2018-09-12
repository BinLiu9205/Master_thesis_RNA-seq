whole_blood<-cbind(summerized$`healthy1,whole_blood`,summerized$`healthy2,whole_blood`,
                   summerized$`healthy3,whole_blood`,summerized$`healthy4,whole_blood`)
Neutrophil<-cbind(summerized$`healthy1,Neutrophil`,summerized$`healthy2,Neutrophil`,
                  summerized$`healthy3,Neutrophil`,summerized$`healthy4,Neutrophil`)
colnames(whole_blood)<-c("healthy_1_whole","healthy_2_whole","healthy_3_whole","healthy_4_whole")
rownames(whole_blood)<-rna_seq$genenames
colnames(Neutrophil)<-c("healthy_1_Neutro","healthy_2_Neutro","healthy_3_Neutro","healthy_4_Neutro")
rownames(Neutrophil)<-rna_seq$genenames

whole_blood_vs_Neutrophil<-cbind(whole_blood,Neutrophil)
library(edgeR)
group<- factor(c(rep("whole_blood",4),rep("Neutrophil",4)))
design<-model.matrix(~0 + group)

y<-DGEList(whole_blood_vs_Neutrophil,group = group)
row.names(y$counts)<-rna_seq[,1]
y <- calcNormFactors(y)
y <- estimateDisp(y, design, robust=TRUE)
fit <- glmQLFit(y, design, robust=TRUE)
#y$common.dispersion
qlf <- glmQLFTest(fit)
whole_blood_Neutrophil<- qlf$table
p.adj<-p.adjust(whole_blood_Neutrophil$PValue,method = "fdr")
whole_blood_Neutrophil$p.adj<-p.adj
