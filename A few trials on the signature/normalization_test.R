## To test whether a DNA/ RNA expressional dataset is normalized

# find the column you are interested in ,it could be any column involving such information
View(summerized)
#have a look at the histogram
hist(rna_seq[,3])
hist(log10(1+rna_seq[,3]))
hist(log10(rna_seq[rna_seq[,3]>0,3]))

View(rna_seq)
hist(apply(summerized,2,mean))
ncol(summerized)
hist(apply(summerized,2,median))
hist(apply(summerized,2,min))
(apply(summerized,2,min))
(apply(summerized,2,median))
nrow(summerized)
hist(apply(summerized,2,sum))
hist(apply(summerized[,names(summerized=="Monocytes")],2,sum))
hist(apply(summerized[,names(summerized)=="Monocytes"],2,sum))
(apply(summerized[,names(summerized)=="Monocytes"],2,sum))
(apply(summerized[,names(summerized)=="Bcells"],2,sum))
(apply(summerized[,names(summerized)=="NK"],2,sum))
(apply(summerized[,names(summerized)=="whole_blood_vs_celltypes.R"],2,sum))
(apply(summerized[,names(summerized)=="whole_blood"],2,sum))
(apply(summerized[,names(summerized)=="whole_blood"],2,var))
apply(summerized[,names(summerized)=="whole_blood"]/,2,var) / apply(summerized[,names(summerized)=="whole_blood"]/,2,sum)
apply(summerized[,names(summerized)=="whole_blood"],2,var) / apply(summerized[,names(summerized)=="whole_blood"]/,2,sum)
apply(summerized[,names(summerized)=="whole_blood"],2,var) / apply(summerized[,names(summerized)=="whole_blood"],2,sum)
1:4 / 1:4
apply(summerized[,names(summerized)=="whole_blood"],1,var) / apply(summerized[,names(summerized)=="whole_blood"],2,sum)
variance_raw <- apply(summerized[,names(summerized)=="whole_blood"],1,var)
head(variance_raw)
sums
sums <- apply(summerized[,names(summerized)=="whole_blood"],2,sum)
cv
cv <- function(x) {sd(x) / mean(x)}
variance_raw <- apply(summerized[,names(summerized)=="whole_blood"],1,sd)
variance_raw <- apply(summerized[,names(summerized)=="whole_blood"],1,cv)
scaled <- summerized / sums
x
x <- matrix(1:9,3)
x
x / c(1,2,3)
x / t(c(1,2,3))
t(t(x) / c(1,2,3))
x
scaled <-  t(t(summerized[,names(summerized)=="whole_blood"]) / sums)
variance_scaled <- apply(scaled,1,cv)
head(scaled)
head(variance_raw)
head(scaled)
head(variance_raw)
head(variance_scaled)
plot(variance_raw,variance_scaled)
abline(a = 0,b = 1,color="red")
abline(a = 0,b = 1,col="red")
lm(variance_raw~variance_scaled)
summary(lm(variance_raw~variance_scaled))