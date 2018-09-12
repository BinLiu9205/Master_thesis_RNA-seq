a_Tol_1 <-read.delim("activated Tol DC _573366.txt",header = FALSE)
a_Tol_2 <-read.delim("activated Tol DC _573367.txt",header = FALSE)
a_Tol_3 <-read.delim("activated Tol DC _573368.txt",header = FALSE)
Tol_1<-read.delim("Tolerogenic_1350464.txt",header = FALSE)
Tol_2<-read.delim("Tolerogenic_1350466.txt",header = FALSE)
Tol_3<-read.delim("tolerogenic DC_573364.txt",header = FALSE)
Tol_4<-read.delim("Tolerogenic DC_573363.txt",header = FALSE)
Tol_5<-read.delim("Tolerogenic DC_573365.txt",header = FALSE)
Tol_6<-read.delim("Tolerogenic DC_1350446.txt",header = FALSE)
Tol_7<-read.delim("Tolerogenic DC_1350450.txt",header = FALSE)
Tol_8<-read.delim("Tolerogenic DC_1350454.txt",header = FALSE)
Tol_9<-read.delim("Tolerogenic DC_1350456.txt",header = FALSE)
Tol_10<-read.delim("Tolerogenic DC_1350463.txt",header = FALSE)
Tol_11<-read.delim("Tolerogenic DC_1350465.txt",header = FALSE)
nTol_1<-read.delim("dex_immature DC_1350444.txt",header = FALSE)
nTol_2<-read.delim("dex_immature DC_1350448.txt",header = FALSE)
nTol_3<-read.delim("dex_immature DC_1350452.txt",header = FALSE)
nTol_4<-read.delim("immature DC_573358.txt",header = FALSE)
nTol_5<-read.delim("immature DC_573359.txt",header = FALSE)
nTol_6<-read.delim("immatureDC_573357.txt",header = FALSE)
nTol_7<-read.delim("immature DC_1350443.txt",header = FALSE)
nTol_8<-read.delim("immature DC_1350447.txt",header = FALSE)
nTol_9<-read.delim("immature DC_1350451.txt",header = FALSE)
nTol_10<-read.delim("immature DC_1350457.txt",header = FALSE)
nTol_11<-read.delim("immature DC_1350458.txt",header = FALSE)
nTol_12<-read.delim("immature DC_1350459.txt",header = FALSE)
nTol_13<-read.delim("mature DC_1350445.txt",header = FALSE)
nTol_14<-read.delim("mature DC_1350449.txt",header = FALSE)
nTol_15<-read.delim("mature DC_1350453.txt",header = FALSE)
nTol_16<-read.delim("mature DC_1350455.txt",header = FALSE)
nTol_17<-read.delim("mature DC_1350460.txt",header = FALSE)
nTol_18<-read.delim("mature DC_1350461.txt",header = FALSE)
nTol_19<-read.delim("mature DC_1350462.txt",header = FALSE)
nTol_20<-read.delim("Maturing DC_573360.txt",header = FALSE)
nTol_21<-read.delim("Maturing DC_573361.txt",header = FALSE)
nTol_22<-read.delim("Maturing DC_573362.txt",header = FALSE)

a_Tol_1<-a_Tol_1[match(toupper(Tol_1[,1]),toupper(a_Tol_1[,1])),]
a_Tol_2<-a_Tol_2[match(toupper(Tol_1[,1]),toupper(a_Tol_2[,1])),]
a_Tol_3<-a_Tol_3[match(toupper(Tol_1[,1]),toupper(a_Tol_3[,1])),]
Tol_3<-Tol_3[match(toupper(Tol_1[,1]),toupper(Tol_3[,1])),]
Tol_4<-Tol_4[match(toupper(Tol_1[,1]),toupper(Tol_4[,1])),]
Tol_5<-Tol_5[match(toupper(Tol_1[,1]),toupper(Tol_5[,1])),]
nTol_4<-nTol_4[match(toupper(Tol_1[,1]),toupper(nTol_4[,1])),]
nTol_5<-nTol_5[match(toupper(Tol_1[,1]),toupper(nTol_5[,1])),]
nTol_6<-nTol_6[match(toupper(Tol_1[,1]),toupper(nTol_6[,1])),]
nTol_20<-nTol_20[match(toupper(Tol_1[,1]),toupper(nTol_20[,1])),]
nTol_21<-nTol_21[match(toupper(Tol_1[,1]),toupper(nTol_21[,1])),]
nTol_22<-nTol_22[match(toupper(Tol_1[,1]),toupper(nTol_22[,1])),]

all_d<-as.data.frame(cbind(Tol_1[,2],Tol_2[,2],Tol_3[,2],Tol_4[,2],
             Tol_5[,2],Tol_6[,2],Tol_7[,2],Tol_8[,2],Tol_9[,2],Tol_10[,2],Tol_11[,2],
             nTol_1[,2],nTol_2[,2], nTol_3[,2],nTol_4[,2],nTol_5[,2],nTol_6[,2],
             nTol_7[,2],nTol_8[,2],nTol_9[,2],nTol_10[,2],nTol_11[,2],nTol_12[,2],
             nTol_13[,2],nTol_14[,2],nTol_15[,2],nTol_16[,2],nTol_17[,2],nTol_18[,2],
             nTol_19[,2],nTol_20[,2],nTol_21[,2],nTol_22[,2]))
colnames(all_d)<-c(rep("Tol",11),rep("nTol",22))

test_d_col <- c(sample(1:11,3),sample(12:33,6))
test_d <- all_d[,test_d_col]
train_d<-all_d[,-c(test_d_col)]
train_d<-train_d[,1:length(train_d[1,])]

m1 <- t(train_d)
d2 <- as.data.frame(m1)
d2$class <- c(rep("Tol",8),rep("nTol",16))

library(ggfortify)
pr_comp <- prcomp(d2[,1:54675])
autoplot(pr_comp, data = d2, colour = 'class')


