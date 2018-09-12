library(limma)
library(AnnotationDbi)
library(org.Hs.eg.db)
healthy_status<-summerized[,grep("sepsis",colnames(summerized))]

healthy_NWB<-healthy_status[,-c(grep("whole_blood",colnames(healthy_status)))]

specific<-factor(c("Neutrophil",rep("Non_Neutrophil",5),"Neutrophil",rep("Non_Neutrophil",4),
                   "Neutrophil",rep("Non_Neutrophil",5)))
design<- model.matrix(~0+specific)

y<-DGEList(healthy_NWB,group =specific)



result<-voom(y,design,plot = T)

fit <- lmFit(result, design)
cont.matrix <- makeContrasts(Compare=specificNeutrophil-specificNon_Neutrophil, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)

fit <- eBayes(fit2)
# result<-topTable(fit, coef=ncol(design),number = 50045,sort.by = "p")

result<-topTable(fit, number = 50045,sort.by = "p")
result_sig<- result[result$adj.P.Val<0.05,]
result_sig$symbol<-mapIds(org.Hs.eg.db,row.names(result_sig),column = "SYMBOL","ENSEMBL")
write.table(result_sig,"sepsis_neutrophil_vs_other_components.txt")
