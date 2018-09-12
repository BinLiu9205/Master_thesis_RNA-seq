library(limma)
View(summerized)

grep("whole_blood",colnames(summerized))
NO_WB<-summerized[,-(c(1,8,15,22,29,36,42))]

status<-factor(c(rep("healthy",24) , rep("sepsis",17)))
component<-factor(c("Neutrophil","Monocytes","Bcells","CD4T","CD8T","NK","Neutrophil","Monocytes","Bcells","CD4T","CD8T","NK",
                    "Neutrophil","Monocytes","Bcells","CD4T","CD8T","NK","Neutrophil","Monocytes","Bcells","CD4T","CD8T","NK",
                    "Neutrophil","Monocytes","Bcells","CD4T","CD8T","NK","Neutrophil","Monocytes","Bcells","CD4T","CD8T",
                    "Neutrophil","Monocytes","Bcells","CD4T","CD8T","NK"))
design<-model.matrix(~0+status+component)
y<-DGEList(NO_WB,group = component)
dge<-estimateDisp(y,design)
plotBCV(dge)

result<-voom(y,design,plot = T)

fit <- lmFit(result, design)
fit <- eBayes(fit)
result<-topTable(fit, coef=ncol(design),number = 50045,sort.by = "p")
result_sig<- result[result$adj.P.Val<0.05,]

write.table(result,file="sepsis_vs_healthy_Neutrophil_single_component.txt",quote = FALSE) 
write.table(result_sig,file="sepsis_vs_healthy_Neutrophil_single_component_significant.txt",quote = FALSE)
