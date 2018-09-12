setwd("D:/00 life science master/masterarbeit R/Project_1_DCs_Glc-DAG induced genes dependant on clec4e")
data_all_symbol<-read.csv('DCs_Glc-DAG induced genes dependant on clec4e.csv')$Gene.Name
library("biomaRt")
library("mygene")
data_all<-queryMany(data_all_symbol,scope="symbol",fields="entrezgene",species="mouse")
datasub_all<-data_all$entrezgene[1:90]
datasub<-datasub_all[is.na(datasub_all)==FALSE]
genes_interested<-c(datasub)
require(goProfiles)
gene.BP <- basicProfile (genes_interested, onto="BP", level=3, orgPackage="org.Mm.eg.db") 
printProfiles(gene.BP, percentage=TRUE)
capture.output(printProfiles(gene.BP, percentage=TRUE),file="DCs_Glc-DAG induced genes dependant on clec4e,goProfile,level=3.txt")
