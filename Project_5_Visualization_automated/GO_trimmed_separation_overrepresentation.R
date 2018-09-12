GO_doc<-read.delim("Panther_result.txt", skip = 12, stringsAsFactors = FALSE) #for Overrepresentation result with skip = 12
m <- regexec(" \\(GO",GO_doc[,1])
GO_doc$trimmed <- substr(GO_doc[,1],1,m)
GO_doc$GO_nr. <- c(substring(GO_doc[,1], regexpr("GO:", GO_doc[,1]) ))
GO_doc$GO_nr. <- substr(GO_doc$GO_nr.,1,10)
GO_doc <- GO_doc[GO_doc$GO_nr. != "Unclassifi",]
GO_doc <- GO_doc[,c(9,1:8)]
colnames(GO_doc) <- c("GO_nr","Term","Size","DE_Genes","Expected","Direction","Enrichment","p_value","GO_name")

write.table(GO_doc$GO_nr, file = "result/GO_number_for_visualization.txt",row.names = F, col.names = F, quote = F)
write.table(GO_doc, file = "result/Panther_result_as_table.txt", row.names = F, quote = F, sep = ";") 
