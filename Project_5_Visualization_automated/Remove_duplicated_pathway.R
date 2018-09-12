all_path<-read.delim("result/GO_best_pathway_for_cytoscape.sif",header = FALSE)
path_sub = all_path[! duplicated(all_path),]
write.table(path_sub, file = "result/best_pathway_nodes.sif",quote=FALSE,col.names = F,row.names = F)
