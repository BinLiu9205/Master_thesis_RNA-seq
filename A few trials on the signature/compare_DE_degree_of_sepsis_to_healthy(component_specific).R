H_spe<-read.table("healthy_people_NK_vs_other_components.txt")
sep_spe<-read.table("sepsis_NK_vs_other_components.txt")

H_only<-H_spe[row.names(H_spe)%in%row.names(sep_spe)==FALSE,]
sep_only<-sep_spe[row.names(sep_spe)%in%row.names(H_spe)==FALSE,]
write.table(H_only,file = "NK_specific_genes_involved_only_in_healthy_people.txt")
write.table(sep_only,file = "NK_specific_genes_involved_only_in_sepsis_people.txt")
H_both<-H_spe[row.names(H_spe)%in%row.names(H_only)==FALSE,]
sep_both<-sep_spe[row.names(sep_spe)%in%row.names(sep_only)==FALSE,]
new<-H_both[match(row.names(sep_both),row.names(H_both)),]
mag<-function(x){2**x}
H_FC<-lapply(new$logFC,mag)
sep_FC<-lapply(sep_both$logFC,mag)
FC_substract<-mapply('-',H_FC,sep_FC,SIMPLIFY=TRUE)
both<-as.data.frame(FC_substract,row.names = row.names(new))

write.table(both,file="differential_DE_extend_for_NK,FChealthy-FCsepsis.txt")
