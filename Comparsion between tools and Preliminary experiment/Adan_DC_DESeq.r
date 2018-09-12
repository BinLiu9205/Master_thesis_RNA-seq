library(DESeq2)
library(ggplot2)
library(ggrepel)
#library(splitstackshape)


setwd("C:/Users/deluca/DATA/Cell Cytometry/DESeq FCs")

data  <- read.delim("Expt174-177.txt",as.is=T,skip = 14)
d <- data[,c("N0174.normalized.","N0175.normalized.","N0176.normalized.","N0177.normalized.")]
names(d) <- c("N0174","N0175","N0176","N0177")
d <- apply(d,2,function(x){as.integer(2^x-1)})
row.names(d) <- data$Gene.ID

#  174: WT OVA, 175 WT Untreated, 176: DKO Untr, 177: DKO Treated (IL17 KO)
info <- data.frame(rbind (c("N0174", "WT OVA"),c("N0175","WT Untr"),c("N0176","DKO Untr"),c("N0177","DKO OVA")))
names(info) <- c("Sample","Condition")

dds <- DESeqDataSetFromMatrix(countData = d, colData = info, design = ~ Condition)
dds <- DESeq(dds,betaPrior=FALSE)
res <- results(dds, contrast=c("Condition","WT OVA","WT Untr"),lfcThreshold=.1, altHypothesis="greaterAbs")

dds.prior <- DESeq(dds)
res.prior.WT <- results(dds.prior, contrast=c("Condition","WT OVA","WT Untr"))
res.prior.KO <- results(dds.prior, contrast=c("Condition","DKO OVA","DKO Untr"))


#bluefunc <- colorRampPalette(c("red", "darkblue"))
plot(res$log2FoldChange,data[,3])
points(res$log2FoldChange[res$padj <= 0.1],data[,3][res$padj <= 0.1], col="red")


plot(res.prior$log2FoldChange,data[,3])
points(res.prior$log2FoldChange[res$padj <= 0.1],data[,3][res$padj <= 0.1], col="red")


plot(res$log2FoldChange,log2(res$baseMean))
plot(res.prior$log2FoldChange,log2(res.prior$baseMean))


annos<- read.delim("../../ensembl_mouse_export2.txt",as.is=T)
res.anno <- res
res.anno <- cbind(res.anno,
                  annos[match(row.names(res),annos$Gene.stable.ID),]
)

res.prior.WT.anno <- res.prior.WT
res.prior.WT.anno <- cbind(res.prior.WT.anno,
                           annos[match(row.names(res.prior.WT),annos$Gene.stable.ID),]
)

res.prior.KO.anno <- res.prior.KO
res.prior.KO.anno <- cbind(res.prior.KO.anno,
                           annos[match(row.names(res.prior.KO),annos$Gene.stable.ID),]
)


plotLines <- function(setName,gsString,annoData,d,info){
  gs <- unlist(strsplit(gsString,","))
  gs <- unique(gs)
  geneNames <- toupper(annoData$Gene.name)
  data.subset <- d[geneNames %in% gs,]
  #names(data.subset)[1]<-"Gene"
  melted <- melt(data.subset)
  melted$Category <- info$Condition[match(melted$Var2,info$Sample)]
  melted$GeneName <- annoData$Gene.name[match(melted$Var1,annoData$Gene.stable.ID)]
  
  ggplot(data=melted, aes(Category,value,fill=GeneName, group=GeneName,color=GeneName)) + geom_line() + scale_y_log10() +
    labs(y="log10 expression",title=setName)
  
}


plotBoxes <- function(setName,gsString,res1,res2,info,label1,label2){
  gs <- unlist(strsplit(gsString,","))
  gs <- unique(gs)
  geneNames <- toupper(res1$Gene.name)
  data.subset1 <- res1[geneNames %in% gs,"log2FoldChange",drop=FALSE]
  data.subset2 <- res2[geneNames %in% gs,"log2FoldChange",drop=FALSE]
  #names(data.subset)[1]<-"Gene"
  melted <- rbind(data.frame(data.subset1,Category = label1, 
                             GeneName= res1$Gene.name[match(row.names(data.subset1),res1$Gene.stable.ID)]),
                  data.frame(data.subset2,Category = label2, 
                             GeneName= res2$Gene.name[match(row.names(data.subset2),res2$Gene.stable.ID)]))
  
  t.pval = t.test(data.subset1$log2FoldChange,data.subset2$log2FoldChange,paired = TRUE)$p.value
  w.pval = wilcox.test(data.subset1$log2FoldChange,data.subset2$log2FoldChange,paired = TRUE)$p.value
  ggplot(data=melted, aes(Category,log2FoldChange)) + geom_boxplot() + 
    labs(y="log2 fold change",title=setName,subtitle=paste("T test: ", t.pval,"; Wilcox:",w.pval)) + 
    geom_jitter(aes(color=GeneName),width = 0.1)  + theme(legend.position="none")
  
}

plotScatter <- function(setName,gsString,res1,res2,info,label1,label2){
  gs <- unlist(strsplit(gsString,","))
  gs <- unique(gs)
  geneNames <- toupper(res1$Gene.name)
  data.subset1 <- res1[geneNames %in% gs,"log2FoldChange",drop=FALSE]
  data.subset2 <- res2[geneNames %in% gs,"log2FoldChange",drop=FALSE]
  #names(data.subset)[1]<-"Gene"
  
  wide <- data.frame(cbind(data.subset1,data.subset2))
  names(wide) <- c(label1,label2)
  wide <- cbind(wide,Gene=res1$Gene.name[match(row.names(data.subset1),res1$Gene.stable.ID)])
  
  ggplot(data=wide,aes(get(label1),get(label2))) + geom_point() + geom_text(aes(label=Gene),hjust=0, vjust=0)+
    geom_vline(xintercept =0) + geom_hline(yintercept =0) + xlab(label1) + ylab(label2) # + geom_abline(slope = 1,intercept = 0)
  
}



setName ="GO_BP_MM_RESPONSE_TO_BIOTIC_STIMULUS"
gsString <- "TUSC5,TMEM91,PRRT2,4930479M11RIK,2200002J24RIK,IFITM7,IFITM5,IFITM2,IFITM1,IFITM3,IFITM6,PRRT1,IFITM10,TMEM233,ACR,GM5331,GM7676"
plotLines(setName,gsString,res.prior.anno,d,info)
plotBoxes(setName,gsString,res.prior.WT.anno,res.prior.KO.anno,info,"WT","KO")
plotScatter(setName,gsString,res.prior.WT.anno,res.prior.KO.anno,info,"WT","KO")

setName ="Type II interferon signaling (IFNG)(WP619) from WikiPathways, data from GeneSetDB"
gsString <- "IRF9,CYBB,ICAM1,IRF8,IFIT2,IFNB1,IFNG,IFNGR1,IFNGR2,IL1B,CXCL10,IRF1,IRF2,IRF4,JAK1,JAK2,CIITA,CXCL9,NOS2,PRKCD,EIF2AK2,PSMB9,PTPN11,REG1,SFPI1,STAT1,STAT2,TAP1,SOCS1,SOCS3,ISG15"
plotLines(setName,gsString,res.prior.WT.anno,d,info)
plotBoxes(setName,gsString,res.prior.WT.anno,res.prior.KO.anno,info,"WT","KO")
plotScatter(setName,gsString,res.prior.WT.anno,res.prior.KO.anno,info,"WT","KO")

setName ="GO_BP_MM_DEFENSE_RESPONSE_TO_GRAM-POSITIVE_BACTERIUM"
gsString <- "C5AR1,TLR2,MBL1,MYD88,IL12A,NOD1,PGLYRP2,PGLYRP1,IL27RA,ACP5,MYO1F,RIPK2,MMP7,GBP7,GBP3,LYZ2,GBP2,LYZ1,HCK,FGR,TIRAP,PLD1,CD36,CTSG,TBK1,NOD2,SCD1,LBP,MBL2,IL6RA,P2RX7,IL6,CAMP,TNFSF8,TNF,LTA,PGLYRP4,PGLYRP3,NCF1,GBP9,GBP10,GBP6,CARD9,KLRK1"
plotLines(setName,gsString,res.prior.anno,d,info)
plotBoxes(setName,gsString,res.prior.WT.anno,res.prior.KO.anno,info,"WT","KO")
plotScatter(setName,gsString,res.prior.WT.anno,res.prior.KO.anno,info,"WT","KO")

setName ="GO_CC_MM_FOCAL_ADHESION"
gsString <- "ZYX,TEK,ITGA5,FERMT2,LIMS1,EZR,LIMD1,LIMS2,GAK,SSH2,SDC4,GIT1,LIMA1,ARHGAP26,PDLIM7,ILK,IRF2,TNS4,VASP,PTPRC,CASS4,PARVA,TLN1,CAV1,APBB1IP,RDX,TENC1,MAPK3,AIF1L,KIF22,MAP2K2,ACTN2,TGFB1I1,SIKE1,PALLD,MYH6,RHOU,MYH7,PDPK1,PIP5K1C,BRCA1,ZFYVE21,FLNB,TRPV4,PTK2,TPPP,HCK,ARHGAP4,PDLIM2,C230081A13RIK,SORBS3,ITGA2B,SORBS1,PXN,LPXN,PARVB,PARVG,MAPK1,MSN,TNS1,LPP,LIG4,ARPC2,FHL2,SDC1,HIC2,EBAG9,PRUNE,BCAR1,EPHA2,STARD8,ARHGEF7,ITGB1,MAP2K1,TNS3,PLEC,FBLIM1,PGM5,DAG1,AATF,ADAM17,FES,VCL,BSPRY,NOX4,PAK1,DLC1,SYNE2,ITGB5,EVL,SH3KBP1,RFWD2,ZFP384,KEAP1,ACTN1,TLN2,ARHGAP24,SORBS2,EPB4.1L5,LUC7L3,ENAH,PNMA1,EPHX2,PTK2B,LIMK1,CDH1,ARHGAP31,MDC1,DIXDC1,CIDEC,TRIP6,LASP1,FERMT1,GRB7"
plotLines(setName,gsString,res.prior.anno,d,info)
plotBoxes(setName,gsString,res.prior.WT.anno,res.prior.KO.anno,info,"WT","KO")
plotScatter(setName,gsString,res.prior.WT.anno,res.prior.KO.anno,info,"WT","KO")

setName ="GO_BP_MM_CELL_MOTILITY"
gsString <- "CD34,ITGB4,SRPX2,MYO10,BRK1,ETS1,TGFBR1,MAP2K1,ENG,ADAM17,DST,FSCN1,RAC1"
plotLines(setName,gsString,res.prior.anno,d,info)
plotBoxes(setName,gsString,res.prior.WT.anno,res.prior.KO.anno,info,"WT","KO")
plotScatter(setName,gsString,res.prior.WT.anno,res.prior.KO.anno,info,"WT","KO")


pathways <- read.delim("../../MousePath_GO_TF_Pathway_Metab.txt")
pathways$genesSym <- as.character(pathways$genesSym)
res1 <- res.prior.WT.anno
res2 <- res.prior.KO.anno
label1 <- "WT"
label2 <- "KO"


geneNames <- toupper(res1$Gene.name)
pathResults <- NULL
for (i in 1:nrow(pathways)){
  #setName <- pathways$label[i]
  gs <- unlist(strsplit(pathways$genesSym[i],","))
  if (length(gs) >= 5 & sum(geneNames %in% gs)>= 5 ){
    #cat("len",length(gs))
    data.subset1 <- res1[geneNames %in% gs,"log2FoldChange",drop=FALSE]
    data.subset2 <- res2[geneNames %in% gs,"log2FoldChange",drop=FALSE]
    
    if(sum(!is.na(data.subset1$log2FoldChange)) >= 5) {
      # signifiance in WT being non-zero
      wt.t.pval = t.test(data.subset1$log2FoldChange)$p.value
      wt.w.pval = wilcox.test(data.subset1$log2FoldChange)$p.value
      # significance in difference between WT and KO
      vs.t.pval = t.test(data.subset1$log2FoldChange,data.subset2$log2FoldChange,paired = TRUE)$p.value
      vs.w.pval = wilcox.test(data.subset1$log2FoldChange,data.subset2$log2FoldChange,paired = TRUE)$p.value
      
      pathResults = rbind(pathResults,data.frame(pathways[i,],wt.t.pval,wt.w.pval,vs.t.pval,vs.w.pval))
    }
  }
}

write.table(pathResults,"pathwayResults.txt",quote=FALSE, sep="\t", row.names = FALSE)


setName ="GO_BP_MM_CELL_DIVISION"
gsString <- "UBE2S,ITGB3BP,SPECC1L,PDCD6IP,DYNC1LI1,WAPAL,SPC25,HAUS4,CENPV,PPP1CB,CDK1,SMC4,DYNLT1A,CCNY,CETN1,DYNLT1B,DYNLT1C,PELO,DYNLT1F,RBBP8,CABLES1,MAEA,ZWINT,HAUS3,MAPRE2,ANAPC16,BUB1B,CDC23,CDC25C,ARL8A,CASC5,UBE2C,NEK2,RGS14,CDCA7,INO80,MAD2L1,CHFR,OIP5,2610002M06RIK,PARD6B,LLGL2,CSNK1A1,WEE1,AURKA,PPP1CA,CEP110,CCNA1,CDK3-PS,SKA1,NEDD1,CHEK2,MAD2L2,PDS5B,NEK9,HAUS1,SEH1L,CDCA5,CTS7,PARD6G,SYCP2,TXNL4A,CHMP1A,CABLES2,BIRC5,SPICE1,SPAG5,MASTL,CETN2,CDK20,RCC1,CDK2,SIRT2,NCAPD3,CDC123,CINP,NEK1,SYCP1,PHF13,SYCE2,KIFC1,NSUN2,HAUS2,MCMBP,FBXO5,CDCA2,STAG1,ARPP19,IST1,HAUS7,BUB3,HAUS5,POGZ,MIS12,PPP2R2D,NUP43,PMF1,LATS1,KIF11,CEP63,RCC2,APITD1,SYCE1,ANAPC13,ANAPC11,MIS18A,CEP55,NDE1,CCNE1,TPX2,ANAPC2,CCNF,STRA13,CENPC1,HELLS,PKN2,ARHGEF2,FZR1,CENPO,ZC3HC1,MAPRE1,E4F1,TACC1,GNAI2,SYCE3,KIF18B,HAUS8,CCND1,NEDD9,CD2AP,NDC80,PPP1R1C,ZFP830,KATNB1,CCND3,LIG3,CENPW,SMC2,NCAPG2,PARD3,CETN3,CDC27,CDK5,TSG101,ZWILCH,NEK6,MAU2,RAD21,CDCA8,CCNG1,SMC1A,MIS18BP1,GNAI3,USP37,NUDC,HMGA2,KIF2C,BIRC6,TIPIN,SMC5,ANAPC4,RB1,KLHL22,GNAI1,PAPD5,SGOL2,RHOA,BOD1,CCNE2,DSN1,SENP5,PTTG1,PSRC1,CDC26,LZTS2,ENSA,TERF1,CDK19,LATS2,CKS1B,RUVBL1,NCAPH,FAM83D,4922501C03RIK,SKA3,CKS2,CCNG2,CDC45,STAG2,VRK1,SGOL1,PPP1CC,CDK7,CDC16,SETDB2,LMLN,NSMCE2,CENPH,DYNLT3,CCNB1,ANAPC7,FSD1,CDK4,SKA2,CDK11B,TEX14,ERCC6L,UBE2I,ANAPC5,RNF8,ANAPC1,CDCA3,KIF20B,VPS4B,PLK5,CDK6,KNTC1,PEX1,MCM5,CDC20,BORA,NCAPD2,KIF2A,PARD6A,CLASP1,KIF2B,SMC3,CENPT,RAN,PARD3B,NEK3,PRKCE,CDK14,MAPRE3,NUP37,CCND2,ANAPC10,CCNA2,ACTR8,TGTP1,CDC14A,CCNB2,ARL8B,CCNO,ZW10,SPC24,CDC7,TIMELESS,CDC25B,CKAP5,NUF2,MAD1L1,CDC6"
plotLines(setName,gsString,res.prior.anno,d,info)
plotBoxes(setName,gsString,res.prior.WT.anno,res.prior.KO.anno,info,"WT","KO")
plotScatter(setName,gsString,res.prior.WT.anno,res.prior.KO.anno,info,"WT","KO")


setName ="GO_BP_MM_INNATE_IMMUNE_RESPONSE"
gsString <- "TLR11,TLR2,MIF,CD55,CD1D1,C4BP,PCBP2,MYD88,TRIL,SPON2,TMEM173,TLR8,TLR7,CD14,TLR13,TICAM2,IL23R,IIGP1,CR1L,CSF1R,PGLYRP1,UNC93B1,IL23A,ECSIT,ZBP1,MALT1,POLR3E,POLR3F,AXL,RIPK2,BTK,MAP4K2,DHX58,IL27,SARM1,POLR3H,SP110,VNN1,POLR3K,IL1RAPL1,IL34,IL1RAPL2,LY86,SIGIRR,CNPY3,CFI,PRKD1,HCK,CHID1,RNF135,TOLLIP,FGR,CSF1,IL1R1,IL1RL2,LYN,IL1RL1,F2RL1,IL18R1,C4B,IL18RAP,FADD,TNFAIP8L2,POLR3D,MST1R,LCN2,PIK3CD,IRAK1,IL1RAP,POLR3G,ITCH,AKIRIN2,BCL10,MAP3K5,BST2,NLRX1,C2,TBK1,TLR12,ISG20,JAK3,DDX58,SAMHD1,KLRG1,TLR9,TLR1,TLR6,LBP,SERPING1,CLEC4A2,CLEC4N,CLEC4D,SYK,POLR3B,C1RL,C1RA,CYBB,LY96,GM5077,IRAK4,CD180,POLR3A,TICAM1,MX1,FCGR1,MX2,TRIM25,TNFSF4,DEFB1,CYBA,TLR3,POLR3C,PGLYRP4,CRCP,TLR4,PGLYRP3,IRGM1,MARCO,IFIH1,CLEC5A,MAVS,TBKBP1"
plotLines(setName,gsString,res.prior.anno,d,info)
plotBoxes(setName,gsString,res.prior.WT.anno,res.prior.KO.anno,info,"WT","KO")
plotScatter(setName,gsString,res.prior.WT.anno,res.prior.KO.anno,info,"WT","KO")


setName ="WIKIPATHWAYS_MM_TOLL-LIKE_RECEPTOR_SIGNALING_PATHWAY-WP75"
gsString <- "AKT3,TLR6,TAB1,TIRAP,CHUK,MAP3K8,MAPK14,TICAM1,AKT1,AKT2,TAB2,FOS,PIK3R5,LY96,TBK1,IFNAR1,IFNAR2,IFNB1,TICAM2,IKBKB,IL1B,IL6,IL12A,IL12B,CXCL10,IRAK1,IRF3,IRF5,IRF7,JUN,LBP,CXCL9,MYD88,NFKB1,NFKB2,NFKBIA,IRAK4,TLR7,TLR8,PIK3CA,PIK3CB,PIK3CD,PIK3CG,PIK3R1,PIK3R2,TLR9,TOLLIP,MAPK1,MAPK3,MAPK8,MAPK11,MAPK9,MAPK10,MAPK13,MAP2K1,MAP2K2,MAP2K3,MAP2K6,MAP2K7,RAC1,RELA,MAPK12,CCL4,CCL5,CXCL11,MAP2K4,SPP1,STAT1,MAP3K7,TLR1,TLR2,TLR3,TLR4,TLR5,TNF,TRAF3,TRAF6,CASP8,PIK3R3,IKBKG,RIPK1,FADD,CD14,CD80,CD86,CD40,IKBKE"
plotLines(setName,gsString,res.prior.anno,d,info)
plotBoxes(setName,gsString,res.prior.WT.anno,res.prior.KO.anno,info,"WT","KO")
plotScatter(setName,gsString,res.prior.WT.anno,res.prior.KO.anno,info,"WT","KO")


setName ="GO_CC_MM_MICROTUBULE"
gsString <- "SLC8A2,TCP1,TRIM54,DYNC1LI1,GM3448,PTPN20,GM3417,TCTE3,RASSF5,HAUS4,HOOK1,TEKT3,TUBA1B,KIF5B,DYNLT1A,DYNLT1B,DYNLT1C,DYNLT1F,DYNC1I2,KIF15,HAUS3,MAPRE2,APC,KIF20A,DNAIC2,KIF19A,NEK2,NAV1,RGS14,INO80,STIM1,NDEL1,KIF21B,CAMSAP2,KIF27,EML2,KIF14,KLC3,MARK4,APOE,KATNAL1,CEP110,KIF16B,TUBA3B,MAP3K11,HSPH1,JAKMIP1,TUBB1,KATNAL2,DNALC4,TUBB6,DNAHC3,HAUS1,TUBB5,SNTB2,CNP,MEFV,DNAHC5,DNAHC17,KIF22,DYNC1H1,SIRT2,CCT5,SARM1,VAPA,TUBB3,KIFC1,HAUS2,STARD9,TBCB,TUBG2,CCT8,CCT3,HAUS7,NDRG1,TEKT5,TUBB4B,TUBB2A,TUBB2B,NINL,CAMSAP3,KIFC5B,HAUS5,GAS8,ODF2,EML3,KLC4,MACF1,KIF11,INCENP,HOOK2,TEKT2,RCC2,STMN1,SPAG16,KIF26A,TUBGCP4,NDE1,TPPP,TPX2,CCT2,VPS41,EML5,KIF3C,TBCA,TTLL9,TTLL1,KIF1B,DNM1,RASSF1,TEKT1,KIF3B,ARHGEF2,KIF5C,DNAHC2,DNM1L,EML6,CCDC105,CEP57,MID2,DCTN1,KIF23,DYNLL1,MAPRE1,CRHBP,DISC1,SPAG6,IQGAP2,KEG1,KIF18B,HAUS8,TPGS1,WDR44,KATNB1,KIFC3,MIPOL1,CSPP1,MAPT,KIF1A,DYNLRB1,MAP1LC3A,GTSE1,FEZ1,CCT7,TUBGCP5,KIF21A,GAS2L2,CAPN6,SPAST,KIF13A,DCX,CENPE,NUDC,TTLL8,MAP2K1,HOOK3,KIF2C,SPAG4,EMD,DYNLRB2,RASSF3,KIF4,TUBA8,TUBE1,NICN1,CENPJ,NCOA2,NIN,SHROOM2,KIF7,HDAC6,DYNC1LI2,PRC1,KIF6,IQGAP1,APPBP2,DYNC1I1,LZTS2,SHROOM3,STAU2,PBXIP1,CEP170,FAM110C,TUBD1,DYNLT3,GABARAP,FSD1,EML1,KIF24,TCP11L1,TUBA4A,KIF26B,KIF12,BCL2L11,MID1IP1,DNAIC1,SLC8A1,EML4,KIF5A,KIF20B,KIFC2,DCTN2,SHROOM1,DNAHC8,TTL,TEKT4,CLIP1,MAP1LC3B,DYNLL2,KIF3A,DYNC2LI1,WDR96,TPPP3,NEIL2,DVL1,KIF2A,LRPPRC,GAS2L1,TUBA3A,DNM3,KIF2B,CAMSAP1,KIF13B,SLC8A3,KIF9,DNM2,CCT6A,MAPRE3,CDK5RAP2,KNCN,FKBP4,POLB,KIF17,CCT4,TUBB4A,GABARAPL1,SPRY2,DYNC2H1,TRIP10,TUBAL3,TTLL3,KIF1C"
plotLines(setName,gsString,res.prior.anno,d,info)
plotBoxes(setName,gsString,res.prior.WT.anno,res.prior.KO.anno,info,"WT","KO")
plotScatter(setName,gsString,res.prior.WT.anno,res.prior.KO.anno,info,"WT","KO")


setName ="GO_MF_MM_MICROTUBULE_MOTOR_ACTIVITY"
gsString <- "DNAHC9,KIF5B,KIF15,KIF20A,DNAIC2,KIF19A,DNAHC1,DNHD1,KIF21B,KIF27,KIF14,KLC3,KIF16B,DNAHC3,KLC2,DNAHC5,DNAHC17,KIF22,DYNC1H1,KIFC1,STARD9,KIFC5B,DNAHC6,KLC1,KLC4,KIF11,KIF26A,KIF3C,KIF1B,KIF3B,KIF5C,DNAHC2,KIF23,DNAHC7A,KIF4-PS,KIF18B,KIFC3,KIF1A,DNAHC11,KIF21A,KIF13A,CENPE,KIF2C,KIF4,KIF7,KIF6,DYNC1I1,DNAHC7B,BBS4,KIF24,KIF26B,KIF12,GM1305,KIF5A,KIF20B,KIFC2,DNAHC8,KIF3A,KIF18A,KIF2A,DNAHC10,KIF2B,KIF13B,DNAHC12,KIF9,KIF17,DYNC2H1,KIF1C"
plotLines(setName,gsString,res.prior.anno,d,info)
plotBoxes(setName,gsString,res.prior.WT.anno,res.prior.KO.anno,info,"WT","KO")
plotScatter(setName,gsString,res.prior.WT.anno,res.prior.KO.anno,info,"WT","KO")

setName ="GO_BP_MM_NEURON_MIGRATION"
gsString <- "CCR4,CXCR4,PCNT,DAB1,CCK,CNTN2,CTNNB1,GAD1,GFRA3,NTN1,NAV1,ATOH1,NDNF,SEMA6A,MYH10,NDEL1,APBB1,SPOCK1,PEX13,NTRK2,MET,ASPM,DCLK1,PRRXL1,VAX1,TBX20,AXL,GPM6A,DCC,PEX2,ALKBH1,MARK2,MATN2,NEUROD4,TUBB2B,TLX3,NEUROG2,SRF,PITX2,KATNA1,ARX,RELN,PTK2,GATA3,NDE1,FEZF1,KIRREL3,DISC1,CXCL12,CDK5R1,NR2F1,BAX,NKX2-1,MAPT,TOP2B,CDK5,PEX7,CELSR1,SATB2,FKTN,DCX,NDN,ROBO3,NR2F2,FYN,NTRK3,LHX1,CCKAR,D130043K22RIK,BARHL1,NR4A2,CELSR2,CDK5R2,DCDC2A,CELSR3,GAS6,PEX5,PAFAH1B1,PRKG1,LMX1B,ESR2,APBB2,MNX1,TWIST1,PAX6,PHOX2B,ASCL1,FZD3,CHL1,PSEN1,GJA1,YWHAE,BARHL2,MARK1,ITGA3"
plotLines(setName,gsString,res.prior.anno,d,info)
plotBoxes(setName,gsString,res.prior.WT.anno,res.prior.KO.anno,info,"WT","KO")
plotScatter(setName,gsString,res.prior.WT.anno,res.prior.KO.anno,info,"WT","KO")

setName ="KEGG_MM_APOPTOSIS"
gsString <- "IRAK2,AKT1,AKT2,APAF1,BIRC3,BIRC2,XIAP,ATM,BAD,BAX,BCL2,BCL2L1,BID,CAPN1,CAPN2,CASP12,CASP3,CASP6,CASP7,CASP8,CASP9,CFLAR,CHUK,CSF2RB,CSF2RB2,CYCS,CYCT,DFFA,DFFB,ENDOG,FADD,FAS,FASL,IKBKB,IKBKG,IL1A,IL1B,IL1R1,IRAK1,IL1RAP,IL3,IL3RA,MYD88,NFKB1,NFKBIA,NGF,NTRK1,PIK3CA,PIK3CD,PIK3R1,PIK3R2,PIK3R3,PRKACA,PRKACB,PPP3CA,PPP3CB,PPP3CC,PPP3R1,PPP3R2,PRKAR1A,PRKAR1B,PRKAR2A,PRKAR2B,PRKX,RELA,RIPK1,TNF,TNFRSF10B,TNFRSF1A,TRAF2,TNFSF10,TRP53,AKT3,IRAK4,AIFM1,PIK3CG,PIK3R5,BIRC7,MAP3K14,CHP1,CHP2,TRADD,IRAK3,PIK3CB"
plotLines(setName,gsString,res.prior.anno,d,info)
plotBoxes(setName,gsString,res.prior.WT.anno,res.prior.KO.anno,info,"WT","KO")
plotScatter(setName,gsString,res.prior.WT.anno,res.prior.KO.anno,info,"WT","KO")

setName ="GO_BP_MM_NEUTROPHIL_CHEMOTAXIS"
gsString <- "C5AR1,ITGA9,ITGB2,NCKAP1L,ITGA1,IL17B,EDN2,EDN3,ITGAM,FCGR3,FCER1G,CX3CL1,IFNG,CCL2,CXCR2,CCL3,CXCL3,CXCL1,CXCL2,CSF3R,CKLF,SLC37A4,PRKCA,SYK,CXADR,ITGB2L,AMICA1,IL1B,SPP1,TGFB2"
plotLines(setName,gsString,res.prior.anno,d,info)
plotBoxes(setName,gsString,res.prior.WT.anno,res.prior.KO.anno,info,"WT","KO")
plotScatter(setName,gsString,res.prior.WT.anno,res.prior.KO.anno,info,"WT","KO")

setName ="GO_BP_MM_CYTOKINESIS"
gsString <- "ROCK1,PIK3C3,KIF20A,RACGAP1,CHMP3,AURKB,ANLN,DSTN,CUL3,PLK1,KLHL9,VPS4A,CFL1,BIRC5,KLHL21,BECN1,INCENP,BCL2L1,RALA,CNTROB,CEP55,RAB11FIP4,PSTPIP1,KLHL13,RAB35,CIT,SON,KIF13A,ROCK2,RAB11A,PRC1,ARL3,AHCTF1,DCTN3,ZFYVE26,MYH9,STX2,RALB,RAB11FIP3"
plotLines(setName,gsString,res.prior.anno,d,info)
plotBoxes(setName,gsString,res.prior.WT.anno,res.prior.KO.anno,info,"WT","KO")
plotScatter(setName,gsString,res.prior.WT.anno,res.prior.KO.anno,info,"WT","KO")

setName ="GO_BP_MM_LYMPHOCYTE_CHEMOTAXIS"
gsString <- "ADAM8,CX3CL1,CCL2,CCL5,CCL3,CKLF"
plotLines(setName,gsString,res.prior.anno,d,info)
plotBoxes(setName,gsString,res.prior.WT.anno,res.prior.KO.anno,info,"WT","KO")
plotScatter(setName,gsString,res.prior.WT.anno,res.prior.KO.anno,info,"WT","KO")


###########################################################################
setName ="Antigen Presentation and Co-stimulation"
gsString <- "Cd40,Cd86,Cd74,Icam1,Swap70,Cd80,Tnfsf9,Tnsfrsf9,Nfkbiz,Tnf"
gsString <- toupper(gsString)
plotLines(setName,gsString,res.prior.anno,d,info)
plotBoxes(setName,gsString,res.prior.WT.anno,res.prior.KO.anno,info,"WT","KO")
plotScatter(setName,gsString,res.prior.WT.anno,res.prior.KO.anno,info,"WT","KO")

setName ="Chemokine Receptors"
gsString <- "Ccr7,Ccr8,ccr5,Ccr4,Ccr9,Ccr11,Ccr12,Ccrl2"
gsString <- toupper(gsString)
plotLines(setName,gsString,res.prior.anno,d,info)
plotBoxes(setName,gsString,res.prior.WT.anno,res.prior.KO.anno,info,"WT","KO")
plotScatter(setName,gsString,res.prior.WT.anno,res.prior.KO.anno,info,"WT","KO")

setName ="Cytokines and Chemokines"
gsString <- "Ccl4,Ccl17,Ccl22,Ccl13,Ccl8,Ccl2,Ccl3,Ccl5,Ccl7,Cx3cl1,Ccrl2,Cxcl9,Cxcl10,Il1a,Il1b,Il6,Il12b"
gsString <- toupper(gsString)
plotLines(setName,gsString,res.prior.anno,d,info)
plotBoxes(setName,gsString,res.prior.WT.anno,res.prior.KO.anno,info,"WT","KO")
plotScatter(setName,gsString,res.prior.WT.anno,res.prior.KO.anno,info,"WT","KO")

setName ="DC Lineage Commitment"
gsString <- "Runx2,Bcl11a,Klf8,Zbtb46,Flt3,Flt3l,Kit,Ccr7,Sfp1,Ikaros,Stat3,Bcl6,Irf8,Batf3,Id2,Mtor,Irf4,Irf2"
gsString <- toupper(gsString)
plotLines(setName,gsString,res.prior.anno,d,info)
plotBoxes(setName,gsString,res.prior.WT.anno,res.prior.KO.anno,info,"WT","KO")
plotScatter(setName,gsString,res.prior.WT.anno,res.prior.KO.anno,info,"WT","KO")

setName ="Antigen Presentation"
gsString <- "Cd74,H2-aa,Icam1,Swap70"
gsString <- toupper(gsString)
plotLines(setName,gsString,res.prior.anno,d,info)
plotBoxes(setName,gsString,res.prior.WT.anno,res.prior.KO.anno,info,"WT","KO")
plotScatter(setName,gsString,res.prior.WT.anno,res.prior.KO.anno,info,"WT","KO")


setName ="Associated with Inflammatory DCs"
gsString <- "Plcl1,Sh3bp4,Gal3st2,Gal3st2,Tnfrsf11a,Lad1,Slamf1,Slco5a1,Ogfrl1,Idh1,Mreg,Rasal2,Rabgap1l,Gpr52,Zbtb37,Rcsd1,Herc4,Adora2a,Efna2,9030607L17Rik,Pctk2,Tmcc3,Gm10752,Plekhg1,Stx11,Myb,1110038D17Rik,Nuak1,Plxnc1,Socs2,Tmtc2,Tmem19,Rassf3,Dtx3,Pisd-ps1,Pus10,Gfpt2,Irf1,Zmynd15,Tm4sf5,Eno3,Arl5c,Cacnb1,Ccr7,Fam49a,Socs2,Tmem18,1700047I17Rik1,Daam1,Gm10432,Actn1,Gtf2a1,Ptpn21,Eml5,Itgb8,Serpinb6b,Hivep1,Aof1,Serinc5,Net1,Sox4,Slc22a23,Cdc14b,Gm10021,B930046C15Rik,3830406C13Rik,E130112L23Rik,Zmym2,Clu,Ccdc122,Epsti1,544988,Gm3264,D830030K20Rik,Hn1l,ENSMUSG00000079376,Hn1l,D830030K20Rik,D830030K20Rik,Gcap14,Psme2,Tbc1d4,Dap,Laptm4b,Apol7e,Gtpbp1,Cacnb3,Spats2,Slc4a8,Ankrd33b,Ankrd33b,Ankrd33b,Rnf19a,Arc,Adcy6,Lima1,Atf7ip2,Snn,Ift57,Cblb,Fam60a,Rogdi,Gm9861,4921513D23Rik,Dnm1l,1700021K19Rik,Ktelc1,Cdgap,Dzip3,Samsn1,Zfp160,Zfp677,Tmem8,Abcg1,H2-gs10,H2-Q8,H2-Q6,Mmp25,H2-T22,H2-M2,Foxp4,Rftn1,Pot1b,Stap2,Arhgap28,Cdc42ep3,Gypc,Fam13b,Pla2g16,Iscu,BC016495,Papss2,Pcgf5,Stxbp3b,AW112010,Frmd4a,Gm10851,Il15ra,Nostrin,Ube2l6,Phf21a,Sema6d,Stard7,Insm1,Stk4,Pfkfb3,Ehmt1,Stxbp1,Traf1,Strbp,Gtdc1,Arl5a,Vps39,Rassf2,Fam164a,Rap2b,Pias3,Rnf115,Gbp5,Gbp2,Gyg,Gnb4,E130311K13Rik,Slc33a1,Nudt17,Polr3c,BC037703,Stxbp3a,Vcam1,Synpo2,N28178,Npr2,Nans,Galnt12,Ttc39a,Pik3r3,Rnf19b,Cdkn2b,Eif2c1,Oprd1,Extl1,Spsb1,Mmp23,Cc2d2a,Anxa3,Bmp2k,Iscu,BC023744,BC023744,BC023744,BC023744,Gatsl2,Gatsl2,Fscn1,Sepsecs,BC062109,5830443L24Rik,BC057170,Gbp4,Mpa2l,Ssh1,Mtif3,N4bp2l1,Fam40b,Trim24,Gimap4,Herc5,A930038C07Rik,Asprv1,Vhl,Rasgef1a,Usp18,Clec2d,Crebl2,Aebp2,B630005N14Rik,Dok1,Htra2,Mxd1,Mical3,Eno2,Tapbpl,Fam60a,4933426I21Rik,Sdhaf1,Il4i1,Mex3b,Ccdc90b,Rab30,Pgap2,Il21r,Prr14,6030429G01Rik,Relb,Kcnk6,1700067C01Rik,Sdhaf1,Lrrk1,Fah,Zfand6,St5,Fam53b,Atp11a,Lamp1,Casp3,D130040H23Rik,Hsh2d,2510049I19Rik,Nlrc5,Nlrc5,Gins3,Atmin,Ap3m2,D030016E14Rik,Il15,Tmem123,Zfp809,Scn4b,Dscaml1,Stoml1,4921528I07Rik,Ube1l,C330006D17Rik,Gramd1b,Sik2,Tspan3,C230081A13Rik,Mpi,4922501C03Rik,4930579C12Rik,Rasa2,E030011O05Rik,Mtmr1,Tbc1d8b,Enox2,Arhgef9,Zmat1,Inpp5d,Cd244,5033414K04Rik,Fuca2,Arhgap18,Mical1,Itgb2,Gna15,Havcr2,Gm2a,Plek,Anxa6,Cd68,Rrm2,Plekhg3,Jarid2,Tgfbi,Dapk1,Mctp1,F2rl2,Gmds,Fgd3,Wdfy2,Gng2,Grap2,Parvg,Nr4a1,Sla,5031439G07Rik,Slc38a1,Racgap1,Itgb7,Ciita,Tfrc,Hlcs,Mapk14,Cyp4f16,Myo1f,H2-DMa,H2-DMb2,H2-Ob,Zfp36l2,Lmnb1,Lpxn,Entpd1,Pik3ap1,Cacnb2,Gsn,Arhgap15,Camk1d,Plcb2,Ehd4,Slc35c2,Nfatc2,Il6ra,BC028528,Cd53,Stmn1,Eif4g3,Sigmar1,Stmn1,Ccnb1,Cenpa,Alox5ap,Sh3tc1,Naaa,Unc119b,Lat2,Hk2,Frmd4b,Plbd1,Tgfb1,Rasgrp4,Napsa,Tm6sf1,Pak1,Nucb2,Itgax,Ptpre,Fosb,Ccnb1,E2f8,Fes,Sult1a1,Stk32c,Jund,Man2b1,Adcy7,Cd209a,Dpep2,Slc9a9,Anln,Plekho2,Ccnb2,Cfp"
gsString <- toupper(gsString)
plotLines(setName,gsString,res.prior.anno,d,info)
plotBoxes(setName,gsString,res.prior.WT.anno,res.prior.KO.anno,info,"WT","KO")
plotScatter(setName,gsString,res.prior.WT.anno,res.prior.KO.anno,info,"WT","KO")

plotScatterDens <- function(setName,gsString,res1,res2,info,label1,label2){
  gs <- unlist(strsplit(gsString,","))
  gs <- unique(gs)
  geneNames <- toupper(res1$Gene.name)
  data.subset1 <- res1[geneNames %in% gs,"log2FoldChange",drop=FALSE]
  data.subset2 <- res2[geneNames %in% gs,"log2FoldChange",drop=FALSE]
  #names(data.subset)[1]<-"Gene"
  
  wide <- data.frame(cbind(data.subset1,data.subset2))
  names(wide) <- c(label1,label2)
  wide <- cbind(wide,Gene=res1$Gene.name[match(row.names(data.subset1),res1$Gene.stable.ID)])
  
  ggplot(data=wide,aes(get(label1),get(label2))) + 
    geom_point(shape = 16, size = 5, show.legend = FALSE, alpha = .4)+ 
    geom_vline(xintercept =0) + geom_hline(yintercept =0) + xlab(label1) + ylab(label2)+ # + geom_abline(slope = 1,intercept = 0)
    geom_text(data=wide[wide[,1] < -1 & wide$Gene!="Dapk1" &wide$Gene!="E2f8"& wide$Gene!="Stmn1",] ,aes(label=Gene),hjust=.5, vjust=-1) +
    geom_text(data=wide[wide[,1] > 1.6& wide$Gene!="Gbp2",] ,aes(label=Gene),hjust=0, vjust=-1)+
    geom_text(data=wide[wide[,2] < -1.5&wide[,1] > -1& wide$Gene!="Oprd1",] ,aes(label=Gene),hjust=-.1, vjust=1) +
    geom_text(data=wide[wide[,2] >1 &wide[,1] < 1& wide$Gene!="Oprd1",] ,aes(label=Gene),hjust=-.1, vjust=-1) +
    ylim(c(-3,3)) + xlim(-4.5,4.5)
  
}
plotScatterDens(setName,gsString,res.prior.WT.anno,res.prior.KO.anno,info,"WT","KO")


#source("https://bioconductor.org/biocLite.R")
#biocLite("globaltest")
library(globaltest)
#gt.options(transpose=TRUE)

#es <- ExpressionSet(d, phenoData=annotatedDataFrameFrom(assayData, byrow=FALSE), featureData=annotatedDataFrameFrom(assayData, byrow=TRUE), experimentData=MIAME(), annotation=character(), protocolData=annotatedDataFrameFrom(assayData, byrow=FALSE), ...)
row.names(info) <- info$Sample
phenoData <- new("AnnotatedDataFrame",data=info)
es <- ExpressionSet(d,phenoData)

gt(Condition, es)

ExpressionSet(assayData, phenoData=annotatedDataFrameFrom(assayData, byrow=FALSE), featureData=annotatedDataFrameFrom(assayData, byrow=TRUE), experimentData=MIAME(), annotation=character(), protocolData=annotatedDataFrameFrom(assayData, byrow=FALSE), ...)

#############################################################################
heatGenes <- NULL
setName ="Antigen Presentation and Co-stimulation"
gsString <- "Cd40,Cd86,Cd74,Icam1,Swap70,Cd80,Tnfsf9,Tnsfrsf9,Nfkbiz,Tnf"
gs <- unlist(strsplit(gsString,","))
heatGenes <- rbind(heatGenes,data.frame(Pathway=setName,Gene=gs))
setName ="Chemokine Receptors"
gsString <- "Ccr7,Ccr8,ccr5,Ccr4,Ccr9,Ccr11,Ccr12,Ccrl2"
gs <- unlist(strsplit(gsString,","))
heatGenes <- rbind(heatGenes,data.frame(Pathway=setName,Gene=gs))
setName ="Cytokines and Chemokines"
gsString <- "Ccl4,Ccl17,Ccl22,Ccl13,Ccl8,Ccl2,Ccl3,Ccl5,Ccl7,Cx3cl1,Ccrl2,Cxcl9,Cxcl10,Il1a,Il1b,Il6,Il12b"
gs <- unlist(strsplit(gsString,","))
heatGenes <- rbind(heatGenes,data.frame(Pathway=setName,Gene=gs))
setName ="DC Lineage Commitment"
gsString <- "Runx2,Bcl11a,Klf8,Zbtb46,Flt3,Flt3l,Kit,Ccr7,Sfp1,Ikaros,Stat3,Bcl6,Irf8,Batf3,Id2,Mtor,Irf4,Irf2"
gs <- unlist(strsplit(gsString,","))
heatGenes <- rbind(heatGenes,data.frame(Pathway=setName,Gene=gs))
setName ="Antigen Presentation"
gsString <- "Cd74,H2-aa,Icam1,Swap70"
gs <- unlist(strsplit(gsString,","))
heatGenes <- rbind(heatGenes,data.frame(Pathway=setName,Gene=gs))

setName ="Associated with Inflammatory DCs"
gsString <- "Plcl1,Sh3bp4,Gal3st2,Gal3st2,Tnfrsf11a,Lad1,Slamf1,Slco5a1,Ogfrl1,Idh1,Mreg,Rasal2,Rabgap1l,Gpr52,Zbtb37,Rcsd1,Herc4,Adora2a,Efna2,9030607L17Rik,Pctk2,Tmcc3,Gm10752,Plekhg1,Stx11,Myb,1110038D17Rik,Nuak1,Plxnc1,Socs2,Tmtc2,Tmem19,Rassf3,Dtx3,Pisd-ps1,Pus10,Gfpt2,Irf1,Zmynd15,Tm4sf5,Eno3,Arl5c,Cacnb1,Ccr7,Fam49a,Socs2,Tmem18,1700047I17Rik1,Daam1,Gm10432,Actn1,Gtf2a1,Ptpn21,Eml5,Itgb8,Serpinb6b,Hivep1,Aof1,Serinc5,Net1,Sox4,Slc22a23,Cdc14b,Gm10021,B930046C15Rik,3830406C13Rik,E130112L23Rik,Zmym2,Clu,Ccdc122,Epsti1,544988,Gm3264,D830030K20Rik,Hn1l,ENSMUSG00000079376,Hn1l,D830030K20Rik,D830030K20Rik,Gcap14,Psme2,Tbc1d4,Dap,Laptm4b,Apol7e,Gtpbp1,Cacnb3,Spats2,Slc4a8,Ankrd33b,Ankrd33b,Ankrd33b,Rnf19a,Arc,Adcy6,Lima1,Atf7ip2,Snn,Ift57,Cblb,Fam60a,Rogdi,Gm9861,4921513D23Rik,Dnm1l,1700021K19Rik,Ktelc1,Cdgap,Dzip3,Samsn1,Zfp160,Zfp677,Tmem8,Abcg1,H2-gs10,H2-Q8,H2-Q6,Mmp25,H2-T22,H2-M2,Foxp4,Rftn1,Pot1b,Stap2,Arhgap28,Cdc42ep3,Gypc,Fam13b,Pla2g16,Iscu,BC016495,Papss2,Pcgf5,Stxbp3b,AW112010,Frmd4a,Gm10851,Il15ra,Nostrin,Ube2l6,Phf21a,Sema6d,Stard7,Insm1,Stk4,Pfkfb3,Ehmt1,Stxbp1,Traf1,Strbp,Gtdc1,Arl5a,Vps39,Rassf2,Fam164a,Rap2b,Pias3,Rnf115,Gbp5,Gbp2,Gyg,Gnb4,E130311K13Rik,Slc33a1,Nudt17,Polr3c,BC037703,Stxbp3a,Vcam1,Synpo2,N28178,Npr2,Nans,Galnt12,Ttc39a,Pik3r3,Rnf19b,Cdkn2b,Eif2c1,Oprd1,Extl1,Spsb1,Mmp23,Cc2d2a,Anxa3,Bmp2k,Iscu,BC023744,BC023744,BC023744,BC023744,Gatsl2,Gatsl2,Fscn1,Sepsecs,BC062109,5830443L24Rik,BC057170,Gbp4,Mpa2l,Ssh1,Mtif3,N4bp2l1,Fam40b,Trim24,Gimap4,Herc5,A930038C07Rik,Asprv1,Vhl,Rasgef1a,Usp18,Clec2d,Crebl2,Aebp2,B630005N14Rik,Dok1,Htra2,Mxd1,Mical3,Eno2,Tapbpl,Fam60a,4933426I21Rik,Sdhaf1,Il4i1,Mex3b,Ccdc90b,Rab30,Pgap2,Il21r,Prr14,6030429G01Rik,Relb,Kcnk6,1700067C01Rik,Sdhaf1,Lrrk1,Fah,Zfand6,St5,Fam53b,Atp11a,Lamp1,Casp3,D130040H23Rik,Hsh2d,2510049I19Rik,Nlrc5,Nlrc5,Gins3,Atmin,Ap3m2,D030016E14Rik,Il15,Tmem123,Zfp809,Scn4b,Dscaml1,Stoml1,4921528I07Rik,Ube1l,C330006D17Rik,Gramd1b,Sik2,Tspan3,C230081A13Rik,Mpi,4922501C03Rik,4930579C12Rik,Rasa2,E030011O05Rik,Mtmr1,Tbc1d8b,Enox2,Arhgef9,Zmat1,Inpp5d,Cd244,5033414K04Rik,Fuca2,Arhgap18,Mical1,Itgb2,Gna15,Havcr2,Gm2a,Plek,Anxa6,Cd68,Rrm2,Plekhg3,Jarid2,Tgfbi,Dapk1,Mctp1,F2rl2,Gmds,Fgd3,Wdfy2,Gng2,Grap2,Parvg,Nr4a1,Sla,5031439G07Rik,Slc38a1,Racgap1,Itgb7,Ciita,Tfrc,Hlcs,Mapk14,Cyp4f16,Myo1f,H2-DMa,H2-DMb2,H2-Ob,Zfp36l2,Lmnb1,Lpxn,Entpd1,Pik3ap1,Cacnb2,Gsn,Arhgap15,Camk1d,Plcb2,Ehd4,Slc35c2,Nfatc2,Il6ra,BC028528,Cd53,Stmn1,Eif4g3,Sigmar1,Stmn1,Ccnb1,Cenpa,Alox5ap,Sh3tc1,Naaa,Unc119b,Lat2,Hk2,Frmd4b,Plbd1,Tgfb1,Rasgrp4,Napsa,Tm6sf1,Pak1,Nucb2,Itgax,Ptpre,Fosb,Ccnb1,E2f8,Fes,Sult1a1,Stk32c,Jund,Man2b1,Adcy7,Cd209a,Dpep2,Slc9a9,Anln,Plekho2,Ccnb2,Cfp"
gs <- unlist(strsplit(gsString,","))
# let us filter this list a bit
wt.ova <- d[match(toupper(gs),toupper(res.prior.WT.anno$Gene.name)),"N0174"]
gs <- gs[!is.na(wt.ova) & wt.ova > 2000]
heatGenes <- rbind(heatGenes,data.frame(Pathway=setName,Gene=gs))

heatGenes$Gene <- toupper(heatGenes$Gene)
#gs <- unique(gs)
geneNames <- toupper(res.prior.WT.anno$Gene.name) # genes as found in the results and data tables
data.subset <- d[match(heatGenes$Gene,geneNames),]
row.names(data.subset) <- make.names(heatGenes$Gene,unique = TRUE)
colnames(data.subset) <- info$Condition

library(pheatmap)
library(RColorBrewer)
cr <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(7)
data.transformed <- na.omit(data.subset[,c("WT Untr","WT OVA", "DKO Untr", "DKO OVA")]+1)
data.transformed <- log2(data.transformed / data.transformed[,1])
breaks <- c(-3,-2,-1,-.5,.5,1,2,3)

anno_col <- heatGenes
row.names(anno_col) <- row.names(data.subset)
anno_col <- anno_col[!is.na(data.subset[,1]),1,drop=FALSE]
pheatmap(data.transformed,show_rownames = FALSE,
         color = cr , #c("blue","lightblue","pink","red"), 
         breaks = breaks,
         cluster_cols = FALSE,cluster_rows = FALSE,
         annotation_row = anno_col, fontsize_row = 5,fontsize_col = 8)

