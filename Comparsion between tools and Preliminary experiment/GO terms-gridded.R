library("ggplot2")
library("ggrepel")
GOdata<-read.delim("shorted.GOtest-result.txt")
GOdata<-GOdata[is.na(GOdata$pvalue_WT)==FALSE&is.na(GOdata$Category)==FALSE,]
GOdata<-GOdata[GOdata$Category=="Immune response"|GOdata$Category=="leukocyte migration",]

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

m <- regexec(" \\(GO",as.character(GOdata$GO.term))
GOdata$trimmed <- substr(GOdata$GO.term,1,m)

pp<-plot<-ggplot(GOdata,aes(x=signedLog(pvalue_WT),y=signedLog(pvalue_KO),color=Significance))+ 
  geom_point()+ xlim(-20,20)+
  geom_text(aes(label=trimmed),size=4) 
pp+facet_wrap( ~ GOdata$Category, ncol=2)

categ <- "Immune response"

ggplot(GOdata[GOdata$Category == categ,],
        aes(x=signedLog(pvalue_WT),y=signedLog(pvalue_KO),color=Significance)
       )+  geom_vline(xintercept = 0) + geom_hline(yintercept = 0) +geom_point()+
  geom_text(aes(label=trimmed),size=3,position=position_jitter(width=1,height=1))

categ <- "leukocyte migration"
ggplot(GOdata[GOdata$Category == categ,],
       aes(x=signedLog(pvalue_WT),y=signedLog(pvalue_KO),color=Significance)
)+  geom_vline(xintercept = 0) + geom_hline(yintercept = 0) +geom_point()+
  geom_text(aes(label=trimmed),size=3,position=position_jitter(width=1,height=1))


categ <- "signaling pathway"
ggplot(GOdata[GOdata$Category == categ,],
       aes(x=signedLog(pvalue_WT),y=signedLog(pvalue_KO),color=Significance)
)+  geom_vline(xintercept = 0) + geom_hline(yintercept = 0) +geom_point()+
  geom_text(aes(label=trimmed),size=3,position=position_jitter(width=0,height=.1))


# everything in one plot, colored by GO term category
ggplot(GOdata,
       aes(x=signedLog(pvalue_WT),y=signedLog(pvalue_KO),color=Category,shape=Significance)
)+  geom_vline(xintercept = 0) + geom_hline(yintercept = 0) +geom_point(size=3)+
 # geom_text(aes(label=trimmed),size=2.5,position=position_jitter(width=1,height=1))+
  geom_text_repel(aes(label=trimmed),size=3,position = position_jitter(width=1,height=1))

ggplot(GOdata,
       aes(x=signedLog(pvalue_WT),y=signedLog(pvalue_KO),color=Category)
)+  geom_vline(xintercept = 0) + geom_hline(yintercept = 0) +geom_point(size=4)

