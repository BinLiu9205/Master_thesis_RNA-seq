nodes<- read.table("result/best_pathway_nodes.sif",as.is=T,header=F)
c_nodes <- as.data.frame(nodes[,1])
colnames(c_nodes)<-"GO_number"
p_nodes <- as.data.frame(nodes[,3])
colnames(p_nodes)<-"GO_number"
A_nodes <- rbind(c_nodes,p_nodes)
A_nodes = A_nodes[! duplicated(A_nodes),]
write.table(A_nodes, file = "result/all_nodes_in the pathway.txt", col.names = F, row.names = F, quote = F)

