add_list <- read.table("result/Panther_result_as_table.txt", as.is = T,header=T, sep= ";")
add_list [1,] <- as.character(colnames(add_list))
add_list
library(xtable)
out_table <- xtable(add_list)
library(jsonlite)
table_js<-toJSON(out_table,dataframe = "rows",pretty = TRUE)
write(paste("var obj =",table_js),file = "result/additional_table.js")
