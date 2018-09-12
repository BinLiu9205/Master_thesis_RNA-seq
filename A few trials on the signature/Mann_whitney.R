d<-read.csv("GSE109125_Gene_count_table.csv")
DC_col<-grep("DC",colnames(d))
DC<-d[,DC_col]
non_DC<-d[,-c(DC_col)]
W_score<-matrix(1,nrow = len_gene,ncol = len_DC)
W_single<-NULL

len_gene<-length(d[,1])
len_DC<-length(DC[1,])
len_all<-length(d[1,])


d[,2:len_all]<-lapply(d[,2:len_all],as.numeric)
DC[,1:len_DC]<-lapply(DC[,1:len_DC],as.numeric)

for (i in 1:len_gene) {
  gene_in<-as.numeric(d[i,2:len_all])
  W_single <- lapply (DC[1,], wilcox.test) $ statistic
  for (j in 1:len_DC) {
    W_score[i,j]<-W_single [j]
}

