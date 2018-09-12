file_list <- list.files(pattern = ".txt")
data_set <- lapply(file_list, function(x){read.delim(x,header = F,as.is = T,stringsAsFactors = F)})
info <- data_set[[1]]
data_set <- unlist(data_set)
data_fr <- as.data.frame(data_set[3:14])
data_fr <- data_fr[,grep("V2", colnames(data_fr))]
#data_set <- lapply(file_list, function(x){read.delim(x,header = F,as.is = T,stringsAsFactors = F)})
rownames(data_fr) <- data_set[3][[1]][,1]
colnames(data_fr) <- info[,2]
matrix_PCA <- t(data_fr)
matrix_PCA <- as.data.frame(matrix_PCA)
matrix_PCA$group <- c("Th2","Th2","Th1","Th17","Th1","Th17","Naive","Naive",
                      "iTreg","iTreg","nTreg","nTreg")

########Absolute intensity (not log transformed) ###############
#########PCA analysis - using all genes as classifier############
library(ggfortify)
library(plotly)
PCA_all <- prcomp(matrix_PCA[,1:45101], center = TRUE, scale. = FALSE)
autoplot(PCA_all, data = matrix_PCA, colour = 'group')
ggplot(PCA_all,aes(PCA_all$x[,1],PCA_all$x[,2],color=matrix_PCA$group)) + geom_point()+
  xlab("PC1") + ylab("PC2") +labs(color="group")
plot(PCA_all,xlab="Dimensions")
PoV <- PCA_all$sdev^2/sum(PCA_all$sdev^2)*100
barplot(PoV, xlab= "Dimensions", ylab="Proportion of explained variance (%)")
#dev.off()
##########PCA analysis - top 500 genes as classifier ##############
top500 <- order(apply(matrix_PCA[,1:45101],2,var),decreasing = TRUE)[1:500]
matrix_top <- matrix_PCA[,top500]
PCA_top500<- prcomp(matrix_top, center = TRUE, scale. = TRUE)
ggplot(NULL,aes(PCA_top500$x[,1],PCA_top500$x[,2],color=matrix_PCA$group)) + geom_point()+
  xlab("PC1") + ylab("PC2") +labs(color="celltypes")
plot(PCA_top500)
PoV_500 <- PCA_top500$sdev^2/sum(PCA_top500$sdev^2)*100
barplot(PoV_500, xlab= "Dimensions", ylab="Proportion of explained variance (%)")
############################################################




Th1 <- data_fr[ ,grep("Th1-",colnames(data_fr),fixed = T)]
Th1 $ mean <- rowMeans(Th1)
Th2 <- data_fr[ ,grep("Th2-",colnames(data_fr),fixed = T)]
Th2 $ mean <- rowMeans(Th2)
Th0 <- data_fr[ ,grep("Naive-",colnames(data_fr),fixed = T)]
Th0$ mean <- rowMeans(Th0)
dir_Th1 <- as.data.frame(Th1$mean / Th0$mean)
dir_Th1 $ direction <- "T"
dir_Th1 $ direction [ dir_Th1$`Th1$mean/Th0$mean` < 1] <- "F"
dir_Th2 <- as.data.frame(Th2$mean / Th0$mean)
dir_Th2 $ direction <- "T"
dir_Th2 $ direction [ dir_Th2$`Th2$mean/Th0$mean` < 1] <- "F"
sum(dir_Th2$direction== dir_Th1$direction)
library(AnnotationDbi)
library(mouse4302.db)
