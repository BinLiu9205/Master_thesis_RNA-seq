setwd("D:/00 life science master/masterarbeit R/Project_1_DCs_Glc-DAG induced genes dependant on clec4e")
data_all<-read.csv('DCs_Glc-DAG induced genes dependant on clec4e.csv')
datasub<-data_all[1:90,]
genenames<-datasub$Gene.Name

library(ontologyIndex)
data(go)

library(ontologySimilarity)
data(gene_GO_terms)
data(GO_IC)

beach <- gene_GO_terms[genenames]
immune<-go$id[go$name=="immune system process"]
beach_immune<-sapply(beach,function(x)intersection_with_descendants(go,roots=immune,x))
data.frame(check.names=FALSE, `#terms`=sapply(beach, length), `#immune terms`=sapply(beach_immune, length))
gene_interested<-beach_immune[sapply(beach_immune, length)!=0]
number<-length(gene_interested)
for (i in 1:number)
{ print(gene_interested[i])
  print(go$name[gene_interested[[i]]])
}

