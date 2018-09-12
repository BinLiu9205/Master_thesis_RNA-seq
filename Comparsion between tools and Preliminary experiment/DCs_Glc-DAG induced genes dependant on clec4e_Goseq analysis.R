library("goseq")
library("geneLenDataBase")
library("org.Mm.eg.db")
library("mygene")
library("GenomicFeatures")
library("biomaRt")
library("ensembldb")
database<-read.table("Mouse Ensembl IDs and Transcript Lengths.txt",header=FALSE,fill=TRUE)
geneset<-read.csv("DCs_Glc-DAG induced genes dependant on clec4e.csv")
#get information about all the genes in the array, save the name of the genes as a character vector
AllGenes<-geneset$Gene.NamegeneLenDataBase
DifGenes<-AllGenes[1:90]
AllGenes.vector<-c(t(as.character(AllGenes)))
DifGenes.vector<-c(t(as.character(DifGenes)))
#find ensemnr of the genes
Ensemb.AllGenes<-queryMany(AllGenes.vector,scope="symbol",fields="ensembl.gene",species="mouse")
Ensemb.DifGenes<-queryMany(DifGenes.vector,scopes="symbol",fields="ensembl.gene",species="mouse")
Ensnr.AllGenes<-Ensemb.AllGenes$ensembl.gene
Ensnr.DifGenes<-Ensemb.DifGenes$ensembl.gene

#gene.vector=as.integer(AllGenes.vector%in%DifGenes.vector))
#produce a txdb database for mouse, automatically if getting access to UCSC mm10
txsByGene=transcriptsBy(TxDb.Mmusculus.UCSC.mm10.ensGene,"gene")
alllengthData=median(width(txsByGene))
txsByGene.frame<-data.frame(txsByGene)
#produce a vector of integer 0 and 1 in order to see whether the gene wanted was expressed
gene.vector=as.integer(Ensnr.AllGenes[is.na(Ensnr.AllGenes)==FALSE]%in%Ensnr.DifGenes[is.na(Ensnr.DifGenes)==FALSE])
#gene.vector=as.integer(AllGenes.vector%in%DifGenes.vector)
alllengthData<-as.data.frame(alllengthData)
#use match to find out all the length data needed for the test
length<-alllengthData[match(Ensnr.AllGenes[is.na(Ensnr.AllGenes)==FALSE],row.names(alllengthData)),]
#Ensnr.AllGenes_sub<-Ensnr.AllGenes[Ensnr.AllGenes%in%rownames(alllengthData)]
#produce a pwf to show how good was the genes fittedgeneLenDataBase"
pwf=nullp(gene.vector,bias.data=length,'ensGene',plot.fit=TRUE)
#produce a dataset of GO terms and add it to the gosgeneLenDataBase"eq function
#mmouse = useDataset(dataset="mmusculus_gene_ensembl",mart=useMart ("ensembl"))
#GOmap = getBM (filters = "ensembl_gene_id", attributes = c("ensembl_gene_id", "go_id"),values= Ensnr.AllGenes,mart = mmouse)
#GOmap<-getgo(Ensnr.AllGenes,'mm10','ensGene',fetch.cats="GO:BP")
#GO.wall<-goseq(pwf,gene2cat = GOmap,'ensGene')
pvals <- goseq(pwf,'mm10','ensGene',test.cats=c("GO:BP"))
#test_pvals<-goseq(pwf,'mm10','ensGene',test.cats=c("GO:BP"),method="Hypergeometric")
#pvals_test<-goseq(pwf,'mm10','knownGene',use_genes_without_cat=TRUE)
enriched.GO=pvals$category[p.adjust(pvals$over_represented_pvalue,method="BH")<.05]
#enriched.GO=pvals$category[p.adjust(pvals_test$over_represented_pvalue,method="BH")<.05]
head(enriched.GO)
#write out the terms interested
library(GO.db)
capture.output(for(go in enriched.GO[1:length(enriched.GO)]){
  print(GOTERM[[go]])
  cat("--------------------------------------\n")
   },file="goseq all result, 14.04.txt")

  