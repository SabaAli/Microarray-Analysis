#   Saba Ali
#   DatatSet Record: GDS2190
#   Title: Bipolar disorder: dorsolateral prefrontal cortex
#   Link: https://www.ncbi.nlm.nih.gov/sites/GDSbrowser?acc=GDS2190

#   Summary: Analysis of postmortem dorsolateral prefrontal cortex from 30 adults with bipolar disorder. 
#            Results provide insight into the pathophysiology of the disease.

#   Platform:	GPL96: [HG-U133A] Affymetrix Human Genome U133A Array


#------------------------------------------#
#    Reading in GEO data and annotations   
#------------------------------------------#
source("https://bioconductor.org/biocLite.R")
biocLite("GEOquery")
library(GEOquery)

gds           <- getGEO("GDS2190")
dat.gds       <- Table(gds)
dat           <- dat.gds
rownames(dat) <- dat[,1] # Makes ID_REF the rownames
dat           <- dat[, -c(1, 2)] # Removes ID_REF & IDENTIFIER column

ann.samples   <- Columns(gds)  # Annotation for samples
ann.genes     <- dat.gds[,1:2] # Annotation for genes/probesets

disease.state <- ann.samples[,2] 
disease.state <- gsub("\\sdisorder",'',disease.state, ignore.case =TRUE) # Removes "disorder" string
disease.state <- factor(disease.state)
colnames(dat) <- paste(colnames(dat), sep = "_", disease.state) 

dim(dat) # There are 61 samples, and 22,283 genes

#------------------------------------------#
#           Removing Outliers
#------------------------------------------#
# Visualize Outliers
library(gplots)
dat.cor <- cor(dat)
dat.avg <- apply(dat.cor,1,mean)

# Average Correlation Plot
par(oma=c(3,0.1,0.1,0.1))
plot(c(1,length(dat.avg)),range(dat.avg),type="n",xlab="",ylab="Avg r",main="Avg correlation of \nBiopolar and Control Subjects",axes=F)
points(dat.avg,bg="red",col=1,pch=21,cex=1.25)
axis(1,at=c(1:length(dat.avg)),labels=dimnames(dat)[[2]],las=2,cex.lab=0.4,cex.axis=0.6)
axis(2)
abline(v=seq(0.5,62.5,1),col="grey")

# Removes potential outliers identified from Average Correlation Plot
dat <- dat[,!names(dat) %in% c("GSM123233_control","GSM123195_bipolar","GSM123204_bipolar","GSM123205_bipolar","GSM123208_bipolar")]
disease.state <- disease.state[c(-1,-61,-60,-59,-58)] # removes outliers from disease.state

dim(dat) # There are 56 samples, and 22,283 genes (after removing outliers)

#------------------------------------------#
#           Filtering Genes
#------------------------------------------#
source("https://bioconductor.org/biocLite.R")
biocLite("affy")
library(affy)

log.dat  <- log2(dat)                  # log2 of data
dat.mean <- apply(log.dat,1,mean)	     # calculate mean for each GENE 
dat.sd   <- sqrt(apply(log.dat,1,var)) # calculate st.deviation for each GENE
dat.cv   <- dat.sd/dat.mean	           # calculate CV for each GENE

#min.cv   <- min(dat.cv) # 0.004273656
#max.cv   <- max(dat.cv) # 0.3380883

dat.25.gene <- quantile(dat.cv,0.25) # 25th quartile

# Plot Histogram of CV values
hist(dat.cv, main="Histogram of CV values \n for Probsets", ylim=c(0,35), cex.main=1.5, col="salmon", prob=T)
abline(v=dat.25.gene,lty=2,col="black") 
lines(density(dat.cv))

# Remove genes with bottom 25% of CV values
q.25           <- quantile(dat.cv,0.25) 
filtered.genes <- dat.cv[dat.cv > q.25] 
length(filtered.genes) # 16,712 (from 22,283)
filtered.genes <- as.data.frame(filtered.genes)
filtered.dat   <- dat[rownames(filtered.genes),] # Subset dat by filtered genes

dat <- filtered.dat
dim(filtered.dat) # 16,712 by 56


#------------------------------------------#
#           Feature Selection
#------------------------------------------#
cl <- colnames(dat)
control <- cl[1:30]
bipolar <- cl[31:56]

# Function to calculate Student’s two-sample t-test on all genes at once
# Returns the p-value for the test
t.test.all.genes <- function(x,s1,s2){
  x1 <- x[s1]
  x2 <- x[s2]
  x1 <- as.numeric(x1)
  x2 <- as.numeric(x2)
  t.out <- t.test(x1, x2, alternative = "two.sided", var.equal = T)
  out <- as.numeric(t.out$p.value)
  return(out)
}

dat.log <- log2(dat)
t.test.run <- apply(dat.log, 1, t.test.all.genes, s1 = control, s2 = bipolar)

# Histogram of P-values from Student's T-Test (Without Feature Selection)
hist(t.test.run, 
     xlab= "p-values",
     main="Histogram of P-values using a Student's t-test \n Bipolar Disorder Study",
     col="lightblue",
     cex.main = 0.9)

sum(t.test.run<0.05) # 1140 probesets (# of genes) with p < 0.05 (threshold)
p.05 <- t.test.run[t.test.run<0.05] # retained genes with p<0.05 with associated p-values
hist(p.05, 
     xlab= "p-values",
     main="Histogram of P-values < 0.05 using a Student's t-test \n Bipolar Disorder Study",
     col="lightgreen",
     cex.main = 0.9)

#------------------------------------------#
#    Clustering/Dimensionality Reduction
#------------------------------------------#
# Subset data by genes determined
length(p.05)
p.05 <- as.data.frame(p.05)
selected.genes <- rownames(p.05)
selected.dat <- dat[selected.genes,]
dim(selected.dat) #1140 x 56

# Principal Component Analysis to visualize samples in 2-D space (First 2 components)
dat.pca <- prcomp(t(selected.dat))
PC1 <- dat.pca$x[,1]
PC2 <- dat.pca$x[,2]
length(PC1)
length(PC2)

plot(PC1,PC2)
plot(range(PC1),range(PC2),type = "n", xlab = "PC1", ylab = "PC2", main = "PCA plot of Bipolar Data" )  
points(PC1[disease.state == "control"], PC2[disease.state == "control"], col=1, bg = "black", pch=21, cex=1.5) 
points(PC1[disease.state == "bipolar"], PC2[disease.state == "bipolar"], col=1, bg = "red", pch=21, cex=1.5) 
legend("bottomright", legend = c("control","bipolar"), pch = 19, col = c("black","red"))  

# Calculate and plot the scree plot that corresponds to the PCA 
dat.pca.var <- round(dat.pca$sdev^2 / sum(dat.pca$sdev^2)*100,2)
plot(c(1:length(dat.pca.var)),
     dat.pca.var,
     type="b",xlab="# Components",ylab="% Variance",pch=21,col=1,bg=3,cex=1.5)
title("Scree plot showing % variability \nexplained by each eigenvalue\nBipolar dataset")

#------------------------------------------#
#               Classification
#------------------------------------------#
library(MASS)
linear.projection <- dat.pca$x[,1:2] # pca data

# First 15 samples of control, first 15 samples of bipolar
train.set <- as.data.frame(rbind(linear.projection[1:15,], linear.projection[33:47,])) # First 15 samples of control, first 15 samples of bipolar
test.set  <- as.data.frame(linear.projection[!(rownames(linear.projection) %in% rownames(train.set)),])
dim(train.set) # 30 x 2

# LDA
train.names <- rownames(train.set)
train.names <- factor(gsub('GSM[[:digit:]]+_', '', train.names))

test.names  <- rownames(test.set)
test.names  <- factor(gsub('GSM[[:digit:]]+_', '', test.names))

lda.train   <- lda(train.names~., train.set)
out         <- predict(lda.train, test.set)
table(out$class,test.names) # confusion matrix

plot(out$x,
     xlab="LD1",ylab="LD2",main="Discriminant function for \nBipolar Dataset",
     bg=as.numeric(test.names), pch=21, col=1)
legend("bottomright", c("Control", "Bipolar"), col = "Black", pt.bg = c("Black", "Red"), pch =21, cex=0.6)

#------------------------------------------#
#     DISCRIMINENT GENES IDENTIFICATION
#------------------------------------------#
ann.probeset  <- dat.gds
rownames(ann.probeset) <- ann.probeset[,1] 
ann.probeset           <- ann.probeset[, -c(1)] 

control.m <- apply(selected.dat[,control],1,mean,na.rm=T) 
bipolar.m <- apply(selected.dat[,bipolar],1,mean,na.rm=T) 

fold <- control.m - bipolar.m
fold.sorted <- sort(fold)

top.5.pos <- tail(fold.sorted,5) # 203798_s_at, 221805_at, 221891_x_at, 208687_x_at, 210338_s_at 
top.5.neg <- head(fold.sorted,5) # 212859_x_at, 211066_x_at, 221501_x_at, 207547_s_at, 204538_x_at

pos.id <- rownames(as.data.frame(top.5.pos))
neg.id <- rownames(as.data.frame(top.5.neg))

pos.identifier <- ann.probeset[pos.id,]
neg.identifier <- ann.probeset[neg.id,]

pos.identifier <- as.data.frame(pos.identifier[,1])
neg.identifier <- as.data.frame(neg.identifier[,1])

rownames(pos.identifier) <- pos.id
rownames(neg.identifier) <- neg.id

pos.identifier
neg.identifier
