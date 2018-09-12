#The database to be compared with (The whole transcriptome related)
alldata<-read.delim("Gene Symbol Log FC ([N0177] vs [N0176])-Treatedvs.Untreated.txt")
#The names of the subset interested in (able to offer 2 alternatives, csv. or txt.)
#genes_to_test<-read.delim("Gene List associated with inflammatory DC.txt")
genes_to_test<-read.csv("DC Lineage commitment gene set.csv")
#get the names we are interested in & remove them from the list
allnames<-as.character(alldata[,1])
names_genes_to_test<-as.character(genes_to_test[,1])
removed_allnames<-allnames[allnames%in%names_genes_to_test==FALSE]
#define two datasets, subdata and maindata, correlate the genes with the Foldchange, or maybe p-value
subdata_to_test<-alldata[match(names_genes_to_test,alldata[,1]),]
subdata_to_test<-subdata_to_test[is.na(subdata_to_test[,1])==FALSE,]
maindata_to_test<-alldata[match(removed_allnames,alldata[,1]),]
maindata_to_test<-maindata_to_test[is.na(maindata_to_test[,1])==FALSE,]

#!!here to define whether to find out the absolute number of the logFC/ p-value, or any parameter could be related to the rank
#recommended as FALSE!
absolute <- FALSE

  if (absolute==TRUE){
    abs_maindata_to_test<-as.data.frame(cbind(as.character(maindata_to_test[,1]),abs(maindata_to_test[,2])))
    abs_subdata_to_test<-as.data.frame(cbind(as.character(subdata_to_test[,1]),abs(as.numeric(subdata_to_test[,2]))))
    result<-wilcox.test(abs_subdata_to_test[,2],as.numeric(abs_maindata_to_test[,2]),alternative="less")
  }else result<-wilcox.test(subdata_to_test[,2],maindata_to_test[,2],alternative="two.sided")

#print result  
print(result)
#could be saved alternative as txt
#capture.output(result,file="result of the Mann whitney test related to defined subset.txt")

