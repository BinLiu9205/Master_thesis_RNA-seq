library("ggplot2")
GOdata<-read.delim("for graph.txt")
GOdata<-GOdata[is.na(GOdata$pvalue_WT)==FALSE,]
WT_only<-subset.data.frame(GOdata,abs(GOdata$pvalue_WT)<0.05&abs(GOdata$pvalue_KO)>0.05)
KO_only<-subset.data.frame(GOdata,abs(GOdata$pvalue_WT)>0.05&abs(GOdata$pvalue_KO)<0.05)
both_signic<-subset.data.frame(GOdata,abs(GOdata$pvalue_WT)<0.05&abs(GOdata$pvalue_KO)<0.05)



GOdata$Significance <- "Both"
GOdata$Significance[abs(GOdata$pvalue_WT) <= 0.05 & abs(GOdata$pvalue_KO) >0.05] <- "WT"
GOdata$Significance[abs(GOdata$pvalue_KO) <= 0.05 & abs(GOdata$pvalue_WT) >0.05] <- "KO"


signedLogSingle <- function (x){
  if (x < 0) log10(-1*x )
  else -log10(x)
}

signedLog <- function (x){
  sapply(x,signedLogSingle)
}

m <- regexec(" \\(GO",GOdata$GO.term)
GOdata$trimmed <- substr(GOdata$GO.term,1,m)

ggplot(GOdata,aes(x=signedLog(pvalue_WT),y=signedLog(pvalue_KO),color=Significance))+ 
  geom_point()+ xlim(-20,20)+
  geom_text(aes(label=trimmed),size=4) 

