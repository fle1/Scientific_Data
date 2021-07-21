## %%%%%%%%%%%%%%%%%%%%%%%%%%%
##
## Script name: Metabolic data pre-processing
##
## Purpose of script: Normalize metabolic data and generation of plots for diagnosis and descriptive analysis
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
## Notes: This script takes as input the raw peak areas resulting from metabolite identification (Provide in supplementary Table 1) and normalize and preprocess the data, resulting in a data.frame ready to be analyzed (Supplementary Table 2)
##
## %%%%%%%%%%%%%%%%%%%%%%%%%%%

# 0. Load table & required libraries ====
library("scales")

metabolites<-read.delim("~/ScientificData/Supplementary_table1.csv",sep = "\t",header = TRUE,comment.char = "#",stringsAsFactors = TRUE)

# Separate the metadata and metabolites in two tables
metadata<-metabolites[,1:7]
metabolites<-metabolites[,8:82]

# Identify each sample by the unique sample name
rownames(metadata)<-metadata$Filename
rownames(metabolites)<-rownames(metadata)

# Delete blank and (technical) control samples
metabolites<-metabolites[-which(metadata$Genotype%in%c("Blank","Arabidopsis leaf control")),]
metadata<-metadata[-which(metadata$Genotype%in%c("Blank","Arabidopsis leaf control")),]
metadata<-droplevels(metadata)

# 1. Normalization ====

## 1.2 Normalize by Internal Standard (Ribitol)
metabolites<-metabolites/metabolites$Ribitol
metabolites<-metabolites[,-which(colnames(metabolites)=="Ribitol")] # Delete Ribitol column

## 1.2 Normalize by sample fresh weight
metabolites<-metabolites/metadata$Weight..mg.

# 2. Check (log) distribution & delete artifacts ====

## Boxplot
boxplot(log(t(metabolites)),outline=FALSE,pch=16,cex=0.8, 
        ylab="log(ug/mgFW)", main = "Per sample metabolite abundance",lwd=0.5,xaxt="n",las=2,
        # Color by Genotype, Tissue and Condition
        col=as.numeric(factor(paste(metadata$Genotype,metadata$Tissue,metadata$Condition))))

# Add lines with the first and second quantiles of the whole dataset
abline(h=quantile(log(metabolites),na.rm = TRUE)[c(2,4)],col="tomato",lty=3,lwd=1.5)


## Overall dataset statistics
quantile(log(metabolites),na.rm=TRUE)

## Identify samples to delete: Sample median, either over the dataset 0.75 percentile or below the dataset 0.25 percentile
wrong_samples<-which(apply(log(metabolites),1,median,na.rm=TRUE) < quantile(log(metabolites),na.rm = TRUE)[2] | 
                         apply(log(metabolites),1,median,na.rm=TRUE) > quantile(log(metabolites),na.rm = TRUE)[4])

## Eliminate as well the duplicated sample 15319oA_1, even though the distribution is not bad
wrong_samples<-c(wrong_samples,which(rownames(metabolites)=="15319oA_1"))

# Print wrong samples, their medians and Q1 & Q3
wrong_samples
apply(log(metabolites[wrong_samples,]),1,median,na.rm=TRUE) # Medians
apply(log(metabolites[wrong_samples,]),1,quantile,na.rm=TRUE) # Quantiles

## Delete these samples for further analysis
metabolites<-metabolites[-wrong_samples,]
metadata<-metadata[-wrong_samples,]

# 3. PCA ====

## Define a function to plot PCA results
myPCA.plot<-function (PCA, main, col = "black", cex = 1, cex.text = 1, pch = 4) {
    P <- PCA$rotation
    TOTvar <- cumsum(PCA$sdev^2/sum(PCA$sdev^2))
    ## Variability due to the two first components
    varPC1 <- round(TOTvar[1] * 100)
    varPC2 <- round((TOTvar[2] - TOTvar[1])*100)
    ## Plot values for the two first components
    plot (P[,1], P[,2], 
          col = col, pch=pch,
          main = main, cex=cex,
          xlim = quantile(P[,1],c(0.01,0.99)),
          ylim = quantile(P[,2],c(0.01,0.99)),
          xlab = paste("PC1: ", varPC1, "% expl.var.", sep = ""),
          ylab = paste("PC2: ", varPC2, "% expl.var.", sep = "")
    )
}


## Calculate principal components per sample (scale the data)
PCA<-prcomp(formula = ~., data = data.frame(t(metabolites)), scale. = TRUE,na.action=na.exclude)

## Calculate principal components only with metabolites that have been detected in both, roots and shoots
tissue_exclusive_mets<-which(apply(metabolites,2,function(x){sum(is.na(x))})>100)
PCA2<-prcomp(formula = ~., data = data.frame(t(metabolites[,-tissue_exclusive_mets])), scale. = TRUE,na.action=na.exclude)

## Calculate principal components only with root sample
PCA3<-prcomp(formula = ~., data = data.frame(t(metabolites[which(metadata$Tissue=="root"),1:49])), scale. = TRUE,na.action=na.exclude)


## Plot PCA. Shape by tissue and Color by drought time point
###  Create a specific vector color to denote watered/time0 (green) or drought (brown to yellow)
drought_colors<-paste(metadata$Condition,metadata$Time..day.,sep = "_")
drought_colors[grep("water",drought_colors)]<-"water"
drought_colors<-factor(drought_colors)
levels(drought_colors)
levels(drought_colors)<-c(terrain.colors(10)[c(4,5,6,7,8,9,1)])

### Plot
myPCA.plot(PCA = PCA,main = "PCA of metabolomic samples",
           pch = c(17,16)[as.numeric(metadata$Tissue)],
           cex = 0.75,
           col = as.character(drought_colors))
### Legends
legend("topleft",fill = terrain.colors(10)[c(1,4,9)], legend = c("time 0 / Water series", "1 day drought", "6 days drought"),bty = "n")
legend("left",pch=c(21,24),bg = "white", legend = c("Shoots","Roots"),bty = "n")



### Plot PCA. Shape by genotype and Color by drought time point
myPCA.plot(PCA = PCA3,main = "PCA of metabolomic samples",
           pch = c(17,16,15)[as.numeric(metadata$Genotype)],
           cex = 0.75,
           col = as.character(drought_colors))
### Legends
legend("topleft",fill = terrain.colors(10)[c(1,4,9)], legend = c("time 0 / Water series", "1 day drought", "6 days drought"),bty = "n")
legend("left",pch=c(21,24,22),bg = "white", legend = c("WT","BRL3ox","quad"),bty = "n")


# 4. Plot profiles of known stress markers (metabolites) ====

mymetabolite<-"Raffinose"
mymedians<-tapply(metabolites[,mymetabolite], list(metadata$Genotype,paste(metadata$Tissue,metadata$Condition,metadata$Time..day.,sep = "_")), median,na.rm=TRUE)
mymedians_WT<-mymedians[3,]

mymedians_WT

## Raffinose in shoots
{
### Series drought
plot(x = c(0:6), y = mymedians_WT[c(20,14,15,16,17,18,19)],
     type = "c",pch= "_", col = "#F0C9C0FF", lwd = 2,
     main = paste(mymetabolite,"levels (Shoots)",sep = " "),xlab = "Time (days)",ylab = paste(mymetabolite, "(ug/mgFW)",sep = " "))
### Medians
points(x = c(0:6), y = mymedians_WT[c(20,14,15,16,17,18,19)], pch = "-", cex = 2.5, col = "black")

### Series water
lines(x = c(0:6), y = mymedians_WT[c(20:26)],
      type = "c",pch= "_", col= "#00A600FF",lwd = 2)
points(x = c(0:6), y = mymedians_WT[c(20:26)], pch = "-", cex = 2.5, col = "black")
    
### Add points
points(x = jitter(as.numeric(subset(metadata,Genotype=="Col-0 WT"&Tissue=="shoot")$Time..day)),
       y = metabolites[rownames(subset(metadata,Genotype=="Col-0 WT"&Tissue=="shoot")),mymetabolite],
       pch = 21, cex = 0.75, col= "gray",bg = alpha(c("#F0C9C0FF","#00A600FF")[as.numeric(subset(metadata,Genotype=="Col-0 WT"&Tissue=="shoot")$Condition)], 0.75))

### Legend
legend("topleft",pch=21, col="black",pt.bg=c("#00A600FF","#F0C9C0FF"),legend=c("Water","Drought"),bty="n")
}

## Raffinose in roots
{
### Series drought
plot(x = c(0:6), y = mymedians_WT[c(7,1:6)],
     type = "c",pch= "_", col = "#F0C9C0FF", lwd = 2,
     main = paste(mymetabolite,"levels (Roots)",sep = " "),xlab = "Time (days)",ylab = paste(mymetabolite, "(ug/mgFW)",sep = " "))
### Medians
points(x = c(0:6), y = mymedians_WT[c(7,1:6)], pch = "-", cex = 2.5, col = "black")

### Series water
lines(x = c(0:6), y = mymedians_WT[c(7:13)],
      type = "c",pch= "_", col= "#00A600FF",lwd = 2)
points(x = c(0:6), y = mymedians_WT[c(7:13)], pch = "-", cex = 2.5, col = "black")

### Add points
points(x = jitter(as.numeric(subset(metadata,Genotype=="Col-0 WT"&Tissue=="root")$Time..day)),
       y = metabolites[rownames(subset(metadata,Genotype=="Col-0 WT"&Tissue=="root")),mymetabolite],
       pch = 21, cex = 0.75, col= "gray",bg = alpha(c("#F0C9C0FF","#00A600FF")[as.numeric(subset(metadata,Genotype=="Col-0 WT"&Tissue=="shoot")$Condition)], 0.75))

### Legend
legend("topleft",pch=21, col="black",pt.bg=c("#00A600FF","#F0C9C0FF"),legend=c("Water","Drought"),bty="n")
}

# 5. Save processed data ====

# Merge metadata and metabolite data frames, and save it as a tab-delimited table
write.table(x = data.frame(metadata,metabolites),file = "Processed_metabolite_table.txt",sep = "\t",col.names = TRUE,row.names = FALSE,quote = FALSE)
# The resulting file should be similar to Supplementary Table 2