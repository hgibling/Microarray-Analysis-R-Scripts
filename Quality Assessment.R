##################################################
##################################################

# Start by Analyzing all 20 arrays from one tissue

TISSUE <- "put type of tissue here"
#for example: TISSUE <- "Liver"


# Have the annotation file (Annotations for Rat Gene 21st.csv) on your desktop

# Have just the Cel files in one folder on your desktop

cel.folder.name <- "put folder name here"
#for example: cel.folder.name <- "Liver Cel Files"
#must type the name exactly as it appears (capitalization, spaces, etc)

##################################################
##################################################


dir.create(paste("~/Desktop", paste(TISSUE, "Microarray Analysis", sep=" "), sep="/"))
main.directory <- paste("~/Desktop", paste(TISSUE, "Microarray Analysis", sep=" "), sep="/")
dir.create(paste(main.directory, "All 20 Arrays", sep="/"))
subdir.all <- paste(main.directory, "All 20 Arrays/", sep="/")

library(oligo)

setwd(paste("~/Desktop", cel.folder.name, sep="/"))

tissue.list <- list.files()
raw.tissue <- read.celfiles(tissue.list)


##### Name Files #####

LC <- paste("LC", 1:5, sep="")
ROS <- paste("ROS", 1:5, sep="")
RSV <- paste("RSV", 1:5, sep="")
ZC <- paste("ZC", 1:5, sep="")
#NOTE: adjust for tissue with mislabeled samples

condition.names <- c(LC, ROS, RSV, ZC)

sampleNames(raw.tissue) <- c(condition.names)
#applies condition names and replicate numbers to each array


##### Quality Assessment #####

plot.colors <- rainbow(20)

dir.create(paste(subdir.all, "QA Images of Raw Data", sep="/"))
subdir.allQA <- paste(subdir.all, "QA Images of Raw Data/", sep="/")


### Boxplot ###

boxplot(raw.tissue, range=1.5, col=plot.colors, xlab="Array", ylab="Log Probe Intensity", main=paste(TISSUE, "Raw Log Probe Intensity", sep=" "))
quartz.save(paste(subdir.allQA, paste(TISSUE, "Raw Boxplot.pdf", sep=" "), sep=""), type="pdf", width=15, height=7)


### Density Plot ###

hist(raw.tissue, col=plot.colors, lty=1, xlab="Log Intensity", ylab="Density", main=paste(TISSUE, "Raw Density Estimation", sep=" "))
legend("topright", inset=0.01, cex=0.75, c(condition.names), col=plot.colors, lty=1)
quartz.save(paste(subdir.allQA, paste(TISSUE, "Raw Density Estimation Plot.pdf", sep=" "), sep=""), type="pdf", width=10, height=7)


### Principal Component Analysis Plot ###

pca.colors <- c(rep("red", 5), rep("green", 5), rep("cyan", 5), rep("purple", 5))
pca.legend.colors <- c("red", "green", "cyan", "purple")
pca.conditions <- c("LC", "ROS", "RSV", "ZC")
pca.numbers <- c(1:5)

raw.expression.matrix <- exprs(raw.tissue)
transposed.raw.expression.matrix <- t(raw.expression.matrix)

pca.values.raw <- prcomp(transposed.raw.expression.matrix)

normal<-par(mar=c(5.1, 4.1, 4.1, 2.1), xpd=F)

par(mar=c(5.1, 4.1, 4.1, 6.1), xpd=T)
#adds space on the side of the graph for the legend

plot(pca.values.raw$x, col=pca.colors, pch=20, main=paste(TISSUE, "Raw PCA Plot", sep=" "))
text(pca.values.raw$x, pos=3, offset=0.2, labels=pca.numbers, cex=0.5)
legend("topright", inset=c(-0.15,0), c(pca.conditions), cex=0.75, col=pca.legend.colors, pch=20)
quartz.save(paste(subdir.allQA, paste(TISSUE, "Raw PCA Plot.pdf", sep=" "), sep=""), type="pdf")
par(normal)
#returns the graph parameters to normal

pca.summary.raw <- summary(pca.values.raw)
proportion.variance.raw <- pca.summary.raw$importance[2:3,1:5]
#gets the values corresponding to the proportion of variance and cumulative variance for the first five principal components

barplot(proportion.variance.raw, beside=T, col=c("black","gray"), main=paste(TISSUE, "Raw Proportion of Variance of Principal Components", sep=" "), xlab="Principal Components", ylab="Percentage")
legend("topleft", inset=0.01, cex=0.75, c("Proportion of Variance", "Cumulative Proportion"), pch=15, col=c("black","gray"))
box()
quartz.save(paste(subdir.allQA, paste(TISSUE, "Raw Proportion of Variance PCA.pdf", sep=" "), sep=""), type="pdf", width=7, height=7)


### Hierarchical Clustering Dendogram ###

transposed <- t(raw.expression.matrix)
distance <- dist(transposed)
sample.clusters <- hclust(distance)

plot(sample.clusters, main=paste(TISSUE, "Raw Hierarchical Cluster Dendogram", sep=" "), xlab="Samples", sub="")
quartz.save(paste(subdir.allQA, paste(TISSUE, "Raw Hierarchical Cluster Dendogram.pdf", sep=" "), sep=""), type="pdf")


### Image Plots ###

dir.create(paste(subdir.allQA, "Image Plots", sep="/"))
subdir.allQA.im <- paste(subdir.allQA, "Image Plots/", sep="/")

for (i in 1:20){
	filename <- paste(TISSUE, i)
	image(raw.tissue[,i], xlab=paste(TISSUE, "Probe Intensity", sep=" "))
	quartz.save(paste(subdir.allQA.im, paste(filename, "Image Plot.png", sep=" "), sep=""), type="png")
}


### Log-Transformed Image Plots ###

dir.create(paste(subdir.allQA, "Log Image Plots", sep="/"))
subdir.allQA.ltim <- paste(subdir.allQA, "Log Image Plots/", sep="/")

for (i in 1:20){
	filename <- paste(TISSUE, i)
	image(raw.tissue[,i], xlab=paste(TISSUE, "Log Probe Intensity", sep=" "), transfo=function(x)x)
	quartz.save(paste(subdir.allQA.ltim, paste(filename, "Log Image Plot.png", sep=" "), sep=""), type="png")
}


### Probe Level Model Fit Quality Assessment ###

raw.tissue.plm <- fitProbeLevelModel(raw.tissue)

dir.create(paste(subdir.allQA, "PLM Plots", sep="/"))
subdir.allQA.plm <- paste(subdir.allQA, "PLM Plots/", sep="/")


# NUSE Plots #

NUSE(raw.tissue.plm, xlab="Array", main=paste(TISSUE, "Normalized Unscaled Standard Errors", sep=" "))
quartz.save(paste(subdir.allQA.plm, paste(TISSUE, "NUSE plot.pdf", sep=" "), sep=""), type="pdf", width=15, height=7)


# RLE Plots #

RLE(raw.tissue.plm, xlab="Array", main=paste(TISSUE, "Relative Log Expression", sep=" "))
quartz.save(paste(subdir.allQA.plm, paste(TISSUE, "RLE plot.pdf", sep=" "), sep=""), type="pdf", width=15, height=7)


# Weights Image Plots #

dir.create(paste(subdir.allQA.plm, "Weights Image Plots", sep="/"))
subdir.allQA.plmW <- paste(subdir.allQA.plm, "Weights Image Plots/", sep="/")

for (i in 1:20){
	filename <- paste(TISSUE, i)
	image(raw.tissue.plm, type="weights", which=i, xlab=paste(TISSUE, "Weights PLM Plot", sep=" "))
	quartz.save(paste(subdir.allQA.plmW, paste(filename, "PLM Weights Image.png", sep=" "), sep=""), type="png")
}


# Residuals Image Plots #

dir.create(paste(subdir.allQA.plm, "Residuals Image Plots", sep="/"))
subdir.allQA.plmR <- paste(subdir.allQA.plm, "Residuals Image Plots/", sep="/")

for (i in 1:20){
	filename <- paste(TISSUE, i)
	image(raw.tissue.plm, type="residuals", which=i, xlab=paste(TISSUE, "Residuals PLM Plot", sep=" "))
	quartz.save(paste(subdir.allQA.plmR, paste(filename, "PLM Residuals Image.png", sep=" "), sep=""), type="png")
}


# Signs of the Residuals Image Plots #

dir.create(paste(subdir.allQA.plm, "Signs of Residuals Image Plots", sep="/"))
subdir.allQA.plmS <- paste(subdir.allQA.plm, "Signs of Residuals Image Plots/", sep="/")

for (i in 1:20){
	filename <- paste(TISSUE, i)
	image(raw.tissue.plm, type="sign.residuals", which=i, xlab=paste(TISSUE, "Signs of the Residuals PLM Plot", sep=" "))
	quartz.save(paste(subdir.allQA.plmS, paste(filename, "PLM Signs Image.png", sep=" "), sep=""), type="png")
}


##### Preprossess the Raw Data #####

normalized.tissue <- rma(raw.tissue)

dir.create(paste(subdir.all, "Preprocessed Data Images", sep="/"))
subdir.all.preproc <- paste(subdir.all, "Preprocessed Data Images/", sep="/")


### Boxplot ###

boxplot(normalized.tissue, range=1.5, col=plot.colors, xlab="Array", ylab="Log Probe Intensity", main=paste(TISSUE, "Preprocessed Log Probe Intensity", sep=" "))
quartz.save(paste(subdir.all.preproc, paste(TISSUE, "Preprossed Boxplot.pdf", sep=" "), sep=""), type="pdf", width=15, height=7)


### Density Plot ###

hist(normalized.tissue, col=plot.colors, lty=1, xlab="Log Intensity", ylab="Density", main=paste(TISSUE, "Preprocessed Density Estimation", sep=" "))
legend("topright", inset=0.01, cex=0.75, c(condition.names), col=plot.colors, lty=1)
quartz.save(paste(subdir.all.preproc, paste(TISSUE, "Preprocessed Density Estimation Plot.pdf", sep=" "), sep=""), type="pdf", width=10, height=7)


### Principal Component Analysis ###

preprocessed.expression.matrix <- exprs(normalized.tissue)
transposed.preprocessed.expression.matrix <- t(preprocessed.expression.matrix)

pca.values.preprocessed <- prcomp(transposed.preprocessed.expression.matrix)

par(mar=c(5.1, 4.1, 4.1, 6.1), xpd=T)
plot(pca.values.preprocessed$x, col=pca.colors, pch=20, main=paste(TISSUE, "Preprocessed PCA Plot", sep=" "))
text(pca.values.preprocessed$x, pos=3, offset=0.2, labels=pca.numbers, cex=0.5)
legend("topright", inset=c(-0.15,0), c(pca.conditions), cex=0.75, col=pca.legend.colors, pch=20)
quartz.save(paste(subdir.all.preproc, paste(TISSUE, "Preprocessed PCA Plot.pdf", sep=" "), sep=""), type="pdf")
par(normal)

pca.summary.preprocessed <- summary(pca.values.preprocessed)
proportion.variance.preprocessed <- pca.summary.preprocessed$importance[2:3,1:5]

barplot(proportion.variance.preprocessed, beside=T, col=c("black","gray"), main=paste(TISSUE, "Preprocessed Proportion of Variance of Principal Components", sep=" "), xlab="Principal Components", ylab="Percentage")
legend("topleft", inset=0.01, cex=0.75, c("Proportion of Variance", "Cumulative Proportion"), pch=15, col=c("black","gray"))
box()
quartz.save(paste(subdir.all.preproc, paste(TISSUE, "Preprocessed Proportion of Variance PCA.pdf", sep=" "), sep=""), type="pdf", width=7, height=7)


### Hierarchical Clustering Dendogram ###

preprocessed.transposed <- t(preprocessed.expression.matrix)
preprocessed.distance <- dist(preprocessed.transposed)
preprocessed.sample.clusters <- hclust(preprocessed.distance)

plot(preprocessed.sample.clusters, main=paste(TISSUE, "Preprocessed Hierarchical Cluster Dendogram", sep=" "), xlab="Samples", sub="")
quartz.save(paste(subdir.all.preproc, paste(TISSUE, "Preprocessed Hierarchical Cluster Dendogram.pdf", sep=" "), sep=""), type="pdf")
