library(edgeR)
d_raw<-read.csv("GSE109125_Gene_count_table.csv")
d_raw[,2:len_all]<-lapply(d_raw[,2:len_all],as.numeric)
d<-edgeR::cpm(d_raw[,2:len_all], normalized.lib.sizes=TRUE, log=FALSE, prior.count=0.25)
d<-as.data.frame(d)

DC_col<-grep("DC",colnames(d))
DC<-d[,c(DC_col)]
len_gene<-length(d[,1])
len_DC<-length(DC_col)
len_all<-length(d[1,])
d[,1:len_all]<-lapply(d[,1:len_all],as.numeric)
DC[,1:len_DC]<-lapply(DC[,1:len_DC],as.numeric)

R_score<-matrix(1,nrow = len_gene , ncol=len_all-1 )
R_single<-NULL

for (i in 1:len_gene){
  R_single <- rank(d[i,1:len_all])
  for (j in 1:(len_all-1)){
    R_score[i,j]<-R_single[j]
  }
}

DC_re<-R_score[,DC_col]
non_DC<-R_score[,-c(DC_col)]
new<-cbind(DC_re,non_DC)
Mann_p<-NULL

for (i in 1:len_gene){
  Mann_p[i] <- wilcox.test(DC_re[i,],non_DC[i,])$p.value
}



potential<-which(Mann_p<0.05)
gene_list<-as.character(d_raw[potential,1]) 


write.table(gene_list,file = "Gene_list_specific_in_dendritic_cells",quote = FALSE)


