healthy<-cbind(summerized$`healthy1,Monocytes`, summerized$`healthy2,Monocytes`,
                          summerized$`healthy3,Monocytes`,summerized$`healthy4,Monocytes`)
sepsis<-cbind(summerized$`sepsis1,Monocytes`,summerized$`sepsis2,Monocytes`,
                 summerized$`sepsis3,Monocytes`)
WB<-cbind(sepsis,healthy)
factor<-c(rep("sepsis",3),rep("healthy",4))
design<-model.matrix(~ 0 + factor)


y<-DGEList(WB,group = factor)
row.names(y$counts)<-rna_seq[,1]
y <- calcNormFactors(y)
y <- estimateDisp(y, design, robust=TRUE)
fit <- glmQLFit(y, design, robust=TRUE)
qlf <- glmQLFTest(fit)
alldata<-qlf$table
alldata$p.adj<-p.adjust(alldata$PValue)

diff<-alldata[alldata$p.adj<0.05,]
write.table(diff,file = "differential expressed genes in sepsis and healthy adults_ neutrophil.txt",quote = FALSE)
