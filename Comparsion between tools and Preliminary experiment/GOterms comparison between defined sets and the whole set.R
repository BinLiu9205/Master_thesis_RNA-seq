GOterms_all<-read.delim("result of Enrichment logFCKO-logFCWT.txt",skip=5)
GO_to_compare<-read.delim("result of Gene List associated with inflammatory DC.txt",skip=7)
same_result<-GOterms_all[GOterms_all[,1]%in%GO_to_compare[,1]==TRUE,]

write.csv(same_result,file = "result.csv",row.names = FALSE)
