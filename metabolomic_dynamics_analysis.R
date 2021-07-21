## %%%%%%%%%%%%%%%%%%%%%%%%%%%
##
## Script name: metabolomics_dynamics_maSigPro
##
## Purpose of script: Analyze differential metabolite dynamics along the drought time course
##
## Author: Fidel Lozano-Elena
##
## Date Created: 2021-07-21
##
## No copyrights on this work (CC0 1.0) [https://creativecommons.org/publicdomain/zero/1.0/]
##
## Email: fidel.lozano@cragenomica.es
##
## %%%%%%%%%%%%%%%%%%%%%%%%%%%
##
## Notes: This script takes as input normalized metabolite abundances (Supplementary Table 2)
##   
##
## %%%%%%%%%%%%%%%%%%%%%%%%%%%


# 0. Load table & required libraries ====
library(maSigPro)

# Read normalized metabolomic data
metabolites<-read.delim("~/ScientificData/Supplementary_table2.csv",sep = "\t",header = TRUE,comment.char = "#",stringsAsFactors = TRUE)

# Separate the metadata and metabolites in two tables
metadata<-metabolites[,1:7]
metabolites<-metabolites[,8:81]
# Identify each sample by the unique sample name
rownames(metadata)<-metadata$Filename
rownames(metabolites)<-rownames(metadata)

# 1. Subset to root tissues and WT vs. BRL3ox under drought ====
metadata_subset<-subset(metadata,Tissue=="root"&Genotype%in%c("Col-0 WT","35S:BRL3-GFP")&c(Condition=="drought"|Time..day.==0))
metabolite_subset<-metabolites[rownames(metadata_subset),]
metabolite_subset<-droplevels(metabolite_subset)

## Delete unidentified (NA) metabolites
which(apply(metabolite_subset,2,function(x){sum(is.na(x))})>40)
metabolite_subset<-metabolite_subset[,-which(apply(metabolite_subset,2,function(x){sum(is.na(x))})>50)]


# 2. Create experimental design matrix & transpose metabolite matrix (required for maSigPro) ====

##Creation of experimental design matrix for maSigPro
experimentaldesign<-data.frame(Time= metadata_subset$Time..day.,
                               # Time point (numerical)
                               Replicate= as.numeric(factor(paste(metadata_subset$Genotype,metadata_subset$Time..day.,sep = "_"))),   
                               #Replicates are not biological replicates, but each of the different conditions
                               WT=as.numeric(factor(metadata_subset$Genotype))-1,
                               # 1 is WT, 0 is not
                               BRL3ox=(as.numeric(factor(metadata_subset$Genotype))-2)*(-1)
                               # 1 is BRL3ox, 0 is not
                               )

rownames(experimentaldesign)<-rownames(metadata_subset)
experimentaldesign
# This matrix might be built manually if user see it clearer

## Transporse the metabolite table to be analyzed
metabolite_subset<-t(metabolite_subset)

# 3. Dynamics analysis - maSigPro ====

## Design
design<-make.design.matrix(experimentaldesign,degree=3)  #Create regression matrix (List), up to polynomic curve of degree 3
design$edesign #Experimental design matrix
design$groups.vector #Regression variables to experimental groups
design$dis #Regression matrix

## Fit design to curves
fit<-p.vector(metabolite_subset, design,Q=0.05, MT.adjust = "BH", min.obs = 20)
summary(fit)
fit$i
fit$FDR
fit$SELEC
fit$p.adjusted

#Find the significant coefficients for each metabolite
tstep<-T.fit(fit,step.method = "backward",alfa = 0.05)
summary(tstep)
tstep$sol #matrix with the p-values with 1-Regression ANOVA, 2-R square of the model adjustment, 3-p-values of the regression coeff. for selected variables
tstep$sig.profiles #expression values for significative metabolites
help(T.fit) #More info

# Get significant genes/metabolites and cluster them
sigs<-get.siggenes(tstep,rsq = 0.55,vars = "groups")
sigs$summary  # Metabolites significatives, well adjusted to the model
sigs$sig.genes$BRL3oxvsWT$g  #Number of metabolites in this group
sigs$sig.genes$BRL3oxvsWT$sig.pvalues

# Plot these clusters profiles
tem<-see.genes(sigs$sig.genes$BRL3oxvsWT, show.fit = T,dis=design$dis,cluster.method = "kmeans",cluster.data = 1,k=4,newX11 = F,legend = F,color.mode = "rainbow")

# Get the metabolites names per cluster
names(tem$cut)<-rownames(sigs$sig.genes$BRL3oxvsWT$coefficients)
tem$cut

# Plot the profile of a particular metabolite or median profile of several of them
PlotGroups(metabolite_subset[c("Raffinose"),], edesign = experimentaldesign,dis = design$dis,groups.vector = design$groups.vector, 
           main = "Raffinose",show.fit = T,alfa = 0.05,legend = F,summary.mode = "median")

PlotGroups(metabolite_subset[c("Proline","Raffinose"),], edesign = experimentaldesign,dis = design$dis,groups.vector = design$groups.vector, 
           main = "Proline & Raffinose",show.fit = T,alfa = 0.05,legend = F,summary.mode = "median")

# 4. Custom plot for maSigPro results ====
yy<-apply(metabolite_subset[c("Proline","Raffinose","Trehalose","Sucrose"),],2,median,na.rm=TRUE)
rm <- matrix(yy, nrow = 1, ncol = length(yy))
rownames(rm) <- c("Median")
colnames(rm) <- rownames(design$dis)
coeffs<-T.fit(rm,design =  design$dis,step.method = "backward",min.obs = 2,alfa = 0.05,nvar.correction = FALSE)$coefficients
coeffs<-as.numeric(coeffs)
colors<-c("steelblue","tomato")

# Plot
plot(x = jitter(design$edesign[,1],factor = 0.5),y=yy,
     type = "p",pch=21,col="black",bg=colors[design$edesign[,4]+1], cex = 0.8,
     main = "Proline, Raffinose, Trehalose & Sucrose cluster", cex.main=0.9,xlab = "days", ylab = "Metabolite accumulation", sub = "Median profile of 4 metabolites",cex.sub=0.8)
# mean line WT
lines(x = unique(design$edesign[,1]),
      y = tapply(yy,INDEX = list(design$edesign[,1],design$edesign[,3]),median,na.rm=TRUE)[,2],
      lty = 1, lwd = 1.5, col = colors[1])
# mean line BRL3ox
lines(x = unique(design$edesign[,1]),
      y = tapply(yy,INDEX = list(design$edesign[,1],design$edesign[,3]),median,na.rm=TRUE)[,1],
      lty = 1, lwd = 1.5, col = colors[2])

# Adjusted polinomic curve line WT
curve(coeffs[1] + coeffs[3] * x + coeffs[5] * x^2 + coeffs[7] * x^3,
      from = 0,add = TRUE, col = colors[1], lwd = 1.5, lty=2)
# Line BRL3 vs WT
curve(coeffs[1] + coeffs[2] + coeffs[3] * x + coeffs[4] * x + coeffs[5] * x^2 + coeffs[6] * x^2 + coeffs[7] * x^3 + coeffs[8] * x ^3,
      from = 0,add = TRUE, col = colors[2], lwd = 1.5, lty=3)

# Legend
legend("topleft", pch = 21, col= "black",pt.bg = colors,legend = c("WT","BRL3ox"),bty="n" )
legend("left",lty = c(1,2),legend = c("Median profile","Best fit"),bty="n", cex = 0.85)

