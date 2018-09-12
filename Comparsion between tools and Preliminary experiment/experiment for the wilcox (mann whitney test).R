allGenes <- 1:100
test <- wilcox.test(allGenes,sample(1:100,10))
test$p.value


null.ps <- NULL
Ws <- NULL
for (i in 1:20000){
  samp <- sample(1:100,5)
  #test <- wilcox.test(samp,allGenes[-samp])
  test <- wilcox.test(allGenes[-samp],samp)
  null.ps <- c(null.ps, test$p.value)
  Ws <- c(Ws,test$statistic)
}
,1:10
hist(null.ps)
hist(Ws)
