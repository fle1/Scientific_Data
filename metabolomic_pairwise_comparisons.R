## %%%%%%%%%%%%%%%%%%%%%%%%%%%
##
## Script name: Metabolomics pairwise comparisons
##
## Purpose of script: Used normalized and curated data to uncover differences in metabolite accumulation at basal conditions
##
## Author: Fidel Lozano-Elena
##
## Date Created: 2021-07-19
##
## No copyrights on this work (CC0 1.0) [https://creativecommons.org/publicdomain/zero/1.0/]
##
## Email: fidel.lozano@cragenomica.es
##
## %%%%%%%%%%%%%%%%%%%%%%%%%%%
##
## Notes: This script takes as input normalized metabolite abundances (Supplementary Table 2)
##
## %%%%%%%%%%%%%%%%%%%%%%%%%%%

# 0. Load table & required libraries ====
library("scales")
## Load custom boxplot functions
source("https://raw.githubusercontent.com/fle1/canolab_scripts/master/metaboliteFC_boxplots.R")

metabolites<-read.delim("~/ScientificData/Supplementary_table2.csv",sep = "\t",header = TRUE,comment.char = "#",stringsAsFactors = TRUE)

## Separate the metadata and metabolites in two tables
metadata<-metabolites[,1:7]
metabolites<-metabolites[,8:81]
## Identify each sample by the unique sample name
rownames(metadata)<-metadata$Filename
rownames(metabolites)<-rownames(metadata)

## Manually transform to nice metabolite names 
#(colnames starting with numbers and/or spaces are transformed when read in R)
colnames(metabolites)
colnames(metabolites)[c(1:5,14,22,28,29,34,46,58,59,66)]<-c("1-6-anhydro-glucose","2-methyl-malate","2-oxoglutarate",
                                                "Hydrobenzoic acid","4-hydroxy-proline","beta-Alanine",
                                                "Fructose-6P","Glucose-1P","Glucose-6P",
                                                "Glycerol-3P","myo-inositol","Ribulose-5P",
                                                "S-methyl-methionine","Sulfuric acid")
colnames(metabolites)

# 1. Identify metabolites differentially accumulated in basal conditions (SHOOTS; t-test) ====
mycomparison<-metadata[which(metadata$Tissue=="shoot" & metadata$Time..day.==0 & metadata$Genotype%in%c("Col-0 WT","35S:BRL3-GFP")),]
mycomparison

## Delete unidentified (NA) metabolites
which(apply(metabolites[rownames(mycomparison),],2,function(x){sum(is.na(x))})>6)
time0<-metabolites[rownames(mycomparison),-which(apply(metabolites[rownames(mycomparison),],2,function(x){sum(is.na(x))})>6)]

## t-test
sigs<-apply(time0,2,function(x){t.test(x[1:5],x[6:10],na.rm=TRUE)$p.value})
sigs<-which(sigs<0.025) # Fairly astringent threshold
sigs

## Boxplot [BRL3ox on top, normalize to WT levels]
par(mar=c(5,7.5,5,2)+0.1) #Bottom, left, top,right
metabolite_FC_boxplots(genotype1_matrix = t(time0[6:10,names(sigs)]),
                       genotype2_matrix = t(time0[1:5,names(sigs)]),
                       control_matrix = t(time0[1:5,names(sigs)]),
                       horizontal = TRUE, color = c("tomato","steelblue","gray65"), 
                       title = "Differentially accumulated metabolites")
legend("bottomright",fill = c("tomato","steelblue","tomato"),legend = c("BRL3ox","Col-0 WT"),bty="n",cex=0.85)



# 2. Ratio shoots/roots ====

## Get an index with sample information (except for tissue)
sample_info<-paste(metadata$Genotype,metadata$Condition,metadata$Time..day.,metadata$Replicate,sep = "_") # Hacer el match con esto.

### Split the metabolites matrix in roots and shoots
shoots<-metabolites[which(metadata$Tissue=="shoot"),]
rownames(shoots)<-sample_info[which(metadata$Tissue=="shoot")]

roots<-metabolites[which(metadata$Tissue=="root"),]
rownames(roots)<-sample_info[which(metadata$Tissue=="root")]

## Get the samples in which measures in roots and shoots are complete. Use only these for the ratio
common_samples<-intersect(rownames(shoots),rownames(roots))

## Get SHOOT/ROOT a matrix with the ratio values
ratio_matrix<-shoots[common_samples,]/roots[common_samples,]

## Compare ratios that are changing between Col-0 WT CTRL & Col-0 WT 6 days drought
rownames(ratio_matrix)
WT_ratios<-ratio_matrix[c(1:5,58:62),]
WT_ratios

### Delete unidentified (NA) metabolites
which(apply(WT_ratios[rownames(WT_ratios),],2,function(x){sum(is.na(x))})>7)
WT_ratios<-WT_ratios[,-which(apply(WT_ratios[rownames(WT_ratios),],2,function(x){sum(is.na(x))})>7)]

sigs<-apply(WT_ratios,2,function(x){t.test(x[1:5],x[6:10])$p.value})
sigs<-which(sigs<0.05)
sigs

## Boxplot [Drought on top, *No relativization to any condition (control_matrix)]
par(mar=c(5,6.5,5,2)+0.1) #Bottom, left, top,right
metabolite_FC_boxplots(genotype1_matrix = t(WT_ratios[6:10,names(sigs)]),
                       genotype2_matrix = t(WT_ratios[1:5,names(sigs)]),
                       control_matrix = t(WT_ratios[1:5,names(sigs)])/t(WT_ratios[1:5,names(sigs)]),
                       horizontal = TRUE,
                       color = c(terrain.colors(3)[2:1],"gray65"), title = "Shoot/Root changes in Col-0 WT",
                       logarithmic = TRUE, ylab = "log(Shoot/Root)")
legend("bottomright",fill = terrain.colors(3)[2:1],legend = c("Drought","Water"),bty="n",cex=0.85)
