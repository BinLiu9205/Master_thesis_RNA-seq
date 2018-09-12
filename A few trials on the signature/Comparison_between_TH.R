filelist= list.files(pattern = ".fpkm_tracking")
datalist = lapply(filelist, function(x)read.delim(x, header = T, as.is = T, stringsAsFactors = F))
datafr = do.call("cbind",datalist)
data_set <- datafr[,grep("FPKM", colnames(datafr))]
data_set <- data_set [, c(1,5,9,13,17,21,25,29,33,37,41,45,49,53,57,61,65,69,73,77)]
colnames (data_set) <- filelist
data_set$gene_id <- datafr[,4]
data_set[,1:21] <- data_set[,c(21,1:20)]
colnames (data_set) <- c("gene_id",filelist)
Th1<- data_set[grep("Th1",colnames(data_set))]
Th17 <- Th1[,-c(1,2,5,6)]
Th1 <- Th1[,c(1,2,5,6)]
Th2 <- data_set[grep("Th2",colnames(data_set))]
Th9 <- data_set[grep("Th9",colnames(data_set))]
Treg <- data_set[grep("Treg",colnames(data_set))]

data_set <- cbind(Th1, Th2, Th9, Th17, Treg)
length_gene <- nrow(data_set)
group <- factor(c(rep("Th1",4),rep("Th2",4),rep("Th9",4),rep("Th17",4),rep("Treg",4)))
expression_level <- as.data.frame(log2(1+data_set))
qqnorm(expression_level[,3])
qqline(expression_level[,3])
colnames(expression_level) <- c("Th1_1A","Th1_1B","Th1_2A","Th1_2B","Th2_1A","Th2_1B","Th2_2A","Th2_2B",
                                "Th9_1A","Th9_1B","Th9_2A","Th9_2B","Th17_1A","Th17_1B","Th17_2A","Th17_2B",
                                "Treg_1A","Treg_1B","Treg_2A","Treg_2B")
res <- cor(expression_level)
correlation_result <- round(res, 4)
library(corrplot)
corrplot(correlation_result, method="circle")
#dev.off()

library(limma)

design <- model.matrix(~ 0 + factor(c(1,1,1,1,2,2,2,2,3,3,3,3,4,4,4,4,5,5,5,5)))
colnames(design) <- c( "Th1","Th2","Th9","Th17","Treg")
fit <- lmFit(expression_level,design =design)
contrast_mat <-makeContrasts(Th2-Th1,Th9-Th1,Th17-Th1,Treg-Th1,levels = design)
fit2 <- contrasts.fit(fit, contrast_mat)
fit <-eBayes(fit2, 0.01)
p.adj<-p.adjust(fit$p.value,method = "fdr")
p.adj<-matrix(as.numeric(p.adj),ncol=4,nrow=23847)

Th1_sig<-as.data.frame(cbind(as.character(datafr$tracking_id),fit$p.value,p.adj))
names(Th1_sig)<-c("Gene_symbol","pvalue_Th2","pvalue_Th9","pvalue_Th17","pvalue_Treg",
                  "padj_Th2","padj_Th9","padj_Th17","padj_Treg")
#Th1_sig <- Th1_sig[(as.numeric(as.character(Th1_sig$padj_Th2)) <0.05 & as.numeric(as.character(Th1_sig$padj_Th9)) <0.05 &as.numeric(as.character(Th1_sig$padj_Th17)) <0.05 &
                     #as.numeric(as.character(Th1_sig$padj_Treg)) <0.05) , ]

Th1_sig <- Th1_sig[(as.numeric(as.character(Th1_sig$padj_Th2)) <0.05 & as.numeric(as.character(Th1_sig$padj_Th9)) <0.05 &as.numeric(as.character(Th1_sig$padj_Th17)) <0.05 ) , ]

library(limma)

design <- model.matrix(~ 0 + factor(c(1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2)))
colnames(design) <- c( "Th1","non_Th1" )
fit <- lmFit(expression_level,design =design)
contrast_mat <-makeContrasts(non_Th1-Th1,levels = design)
fit2 <- contrasts.fit(fit, contrast_mat)
fit <-eBayes(fit2, 0.01)
p.adj<-p.adjust(fit$p.value,method = "fdr")
p.adj<-matrix(as.numeric(p.adj))

Th1_sig<-as.data.frame(cbind(as.character(datafr$tracking_id),fit$p.value,p.adj))
names(Th1_sig)<-c("Gene_symbol","pvalue_non_Th1","padj_non_Th1")
#Th1_sig <- Th1_sig[(as.numeric(as.character(Th1_sig$padj_Th2)) <0.05 & as.numeric(as.character(Th1_sig$padj_Th9)) <0.05 &as.numeric(as.character(Th1_sig$padj_Th17)) <0.05 &
#as.numeric(as.character(Th1_sig$padj_Treg)) <0.05) , ]

Th1_sig <- Th1_sig[as.numeric(as.character(Th1_sig$padj_non_Th1)) <0.05,]
#################################

design <- model.matrix(~ 0 + factor(c(1,1,1,1,2,2,2,2,3,3,3,3,4,4,4,4,5,5,5,5)))
colnames(design) <- c( "Th1","Th2","Th9","Th17","Treg")
fit <- lmFit(expression_level,design =design)
contrast_mat <-makeContrasts(Th1-Th2,Th9-Th2,Th17-Th2,Treg-Th2,levels = design)
fit2 <- contrasts.fit(fit, contrast_mat)
fit <-eBayes(fit2, 0.01)
p.adj<-p.adjust(fit$p.value,method = "fdr")
p.adj<-matrix(as.numeric(p.adj),ncol=4,nrow=23847)

Th2_sig<-as.data.frame(cbind(as.character(datafr$tracking_id),fit$p.value,p.adj))
names(Th2_sig)<-c("Gene_symbol","pvalue_Th1","pvalue_Th9","pvalue_Th17","pvalue_Treg",
                  "padj_Th1","padj_Th9","padj_Th17","padj_Treg")
#Th1_sig <- Th1_sig[(as.numeric(as.character(Th1_sig$padj_Th2)) <0.05 & as.numeric(as.character(Th1_sig$padj_Th9)) <0.05 &as.numeric(as.character(Th1_sig$padj_Th17)) <0.05 &
#as.numeric(as.character(Th1_sig$padj_Treg)) <0.05) , ]

Th2_sig <- Th2_sig[(as.numeric(as.character(Th2_sig$padj_Th1)) <0.05 & as.numeric(as.character(Th2_sig$padj_Th9)) <0.05 &as.numeric(as.character(Th2_sig$padj_Th17)) <0.05 ) , ]

###########################apply this comparison to the Th2 type##############
library(limma)

design <- model.matrix(~ 0 + factor(c(2,2,2,2,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2)))
colnames(design) <- c( "Th2","non_Th2" )
fit <- lmFit(expression_level,design =design)
contrast_mat <-makeContrasts(non_Th2-Th2,levels = design)
fit2 <- contrasts.fit(fit, contrast_mat)
fit <-eBayes(fit2, 0.01)
p.adj<-p.adjust(fit$p.value,method = "fdr")
p.adj<-matrix(as.numeric(p.adj))

Th2_sig<-as.data.frame(cbind(as.character(datafr$tracking_id),fit$p.value,p.adj))
names(Th2_sig)<-c("Gene_symbol","pvalue_non_Th2","padj_non_Th2")
#Th1_sig <- Th1_sig[(as.numeric(as.character(Th1_sig$padj_Th2)) <0.05 & as.numeric(as.character(Th1_sig$padj_Th9)) <0.05 &as.numeric(as.character(Th1_sig$padj_Th17)) <0.05 &
#as.numeric(as.character(Th1_sig$padj_Treg)) <0.05) , ]

Th2_sig <- Th2_sig[as.numeric(as.character(Th2_sig$padj_non_Th2)) <0.05,]

##########Treg as a trial################
library(limma)

design <- model.matrix(~ 0 + factor(c(2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,1,1,1)))
colnames(design) <- c( "Treg","non_Treg" )
fit <- lmFit(expression_level,design =design)
contrast_mat <-makeContrasts(non_Treg-Treg,levels = design)
fit2 <- contrasts.fit(fit, contrast_mat)
fit <-eBayes(fit2, 0.01)
p.adj<-p.adjust(fit$p.value,method = "fdr")
p.adj<-matrix(as.numeric(p.adj))

Treg_sig<-as.data.frame(cbind(as.character(datafr$tracking_id),fit$p.value,p.adj))
names(Treg_sig)<-c("Gene_symbol","pvalue_non_Treg","padj_non_Treg")
#Th1_sig <- Th1_sig[(as.numeric(as.character(Th1_sig$padj_Th2)) <0.05 & as.numeric(as.character(Th1_sig$padj_Th9)) <0.05 &as.numeric(as.character(Th1_sig$padj_Th17)) <0.05 &
#as.numeric(as.character(Th1_sig$padj_Treg)) <0.05) , ]

Treg_sig <- Treg_sig[as.numeric(as.character(Treg_sig$padj_non_Treg)) <0.05,]


