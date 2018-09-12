#trial

trial<-read.table("sepsis_vs_healthy_Bcells_significant.txt")
all_over<-read.table("sepsis_vs_healthy_Neutrophil_single_component_significant.txt")

unmatched_1<-trial[row.names(trial)%in%row.names(all_over)==FALSE,]


write.table(unmatched_1,file = "sepsis_vs_healthy_Bcells_specific_significant.txt")

