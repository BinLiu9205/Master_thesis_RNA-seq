whole_blood<-read.table("sepsis_vs_healthy_whole_blood_significant.txt")
Neutrophil<-read.table("sepsis_vs_healthy_Neutrophil_significant.txt")
Monocytes<-read.table("sepsis_vs_healthy_Monocytes_significant.txt")
Bcells<-read.table("sepsis_vs_healthy_Bcells_significant.txt")
CD4T<-read.table("sepsis_vs_healthy_CD4T_significant.txt")
CD8T<-read.table("sepsis_vs_healthy_CD8T_significant.txt")
NK<-read.table("sepsis_vs_healthy_NK_significant.txt")

unique<-whole_blood[row.names(whole_blood)%in%row.names(Neutrophil)==FALSE&&row.names(whole_blood)%in%row.names(Monocytes)==FALSE&&
                    row.names(whole_blood)%in%row.names(Bcells)==FALSE&&row.names(whole_blood)%in%row.names(CD4T)==FALSE&&
                    row.names(whole_blood)%in%row.names(CD8T)==FALSE&&row.names(whole_blood)%in%row.names(NK)==FALSE,]
write.table(unique,file = "DE genes in whole_blood but not in immune related subsets.txt")

