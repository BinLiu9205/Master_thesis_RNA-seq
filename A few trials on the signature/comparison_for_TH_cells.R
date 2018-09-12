file_list <- list.files(pattern = ".txt")
data_set <- lapply(file_list, function(x){read.delim(x,header = F,as.is = T,stringsAsFactors = F)})
info <- data_set[[1]]
unlist(data_set)
data_fr <- as.data.frame(data_set[3:14])
data_fr <- data_fr[,grep("V2", colnames(data_fr))]
rownames(data_fr) <- data_set[3][[1]][,1]
colnames(data_fr) <- info[,2]
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
dir_Th1$symbol <- mapIds(mouse4302.db, keys = row.names(data_fr), column = "SYMBOL",keytype = "PROBEID")
dir_Th2$symbol <- mapIds(mouse4302.db,keys= row.names(data_fr),column = "SYMBOL",keytype = "PROBEID")
diff_genes <- read.table("Diff_genes_Th1_Th2.txt",header =T) 
sig_Th1 <- dir_Th1[match(diff_genes$Gene_symbol,dir_Th1$symbol),]
sig_Th1 <- sig_Th1[is.na(sig_Th1$symbol)==F,]
sig_Th2 <- dir_Th2[match(diff_genes$Gene_symbol,dir_Th2$symbol),]
sig_Th2 <- sig_Th2[is.na(sig_Th2$symbol)==F,]
sum(sig_Th2$direction!=sig_Th1$direction)
diff_direction <- sig_Th1$symbol[sig_Th2$direction!=sig_Th1$direction]
Th1_diff <- sig_Th1[sig_Th2$direction!=sig_Th1$direction,]
Th2_diff <- sig_Th2[sig_Th2$direction!=sig_Th1$direction,]
diff_direction

#########
