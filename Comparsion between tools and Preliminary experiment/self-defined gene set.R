alldata<-read.delim("logKO-logWT_data.txt")
data_wanted<-read.delim("Gene List associated with inflammatory DC.txt")
data_wanted_logFC<-alldata[match(data_wanted[,1],alldata[,1]),]
data_wanted_logFC<-data_wanted_logFC[is.na(data_wanted_logFC$log.177vs176..log.174vs175.)==FALSE,]
write.table(data_wanted_logFC,file="Gene_related_to_inflammatory_DC_logFC.txt")

