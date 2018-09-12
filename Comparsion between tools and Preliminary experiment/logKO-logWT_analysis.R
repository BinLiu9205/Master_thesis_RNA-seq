logdata<-read.delim("Enrichment logFCKO-logFCWT.txt",skip=5)
Compare<-read.delim("logKO-logWT_hierarchical.txt")
new_logdata<-logdata[match(Compare$hierarchical,logdata$GO.biological.process.complete),]
hierarchical<-new_logdata[is.na(new_logdata$pvalue)==FALSE,]
hierarchical$regulation="upregulated"
hierarchical$regulation[hierarchical$overUnder=="-"]<-"downregulated"
others<-logdata[logdata$GO.biological.process.complete%in%hierarchical$GO.biological.process.complete==FALSE,]
others$regulation="upregulated"
others$regulation[others$overUnder=="-"]<-"downregulated"
write.table(hierarchical,file="GOanalysis_logFC-minus.txt",col.names = FALSE,row.names = FALSE)
write.table(others,file="GOanalysis_logFC-minus.txt",col.names = FALSE,row.names = FALSE,append = TRUE)
