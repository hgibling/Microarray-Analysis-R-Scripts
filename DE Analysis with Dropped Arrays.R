##################################################
##################################################

# If after performing Quality Assessment it is decided some arrays should be dropped: 

# Have the annotation file (Annotations for Rat Gene 21st.csv) on your desktop

arrays.to.be.dropped <- c(0, 0, 0)
#for example: arrays.to.be.dropped<-c(7,13)

# Input remaining arrays for each condition

LC.remaining <- c(0:0)
ROS.remaining <- c(0:0)
RSV.remaining <- c(0:0)
ZC.remaining <- c(0:0)
#for example: ROS.remaining <- c(1, 3:5)
#if array 7 has been dropped, the remaining ROS samples are ROS1, ROS3, ROS4 and ROS5

##################################################
##################################################


dir.create(paste(main.directory, "Data With Dropped Arrays", sep="/"))
subdir.drop <- paste(main.directory, "Data With Dropped Arrays/", sep="/")

raw.tissue.dropped <- raw.tissue[,-c(arrays.to.be.dropped)]

LC.names <- paste("LC", LC.remaining, sep="")
ROS.names <- paste("ROS", ROS.remaining, sep="")
RSV.names <- paste("RSV", RSV.remaining, sep="")
ZC.names <- paste("ZC", ZC.remaining, sep="")

condition.names.dropped <- c(LC.names, ROS.names, RSV.names, ZC.names)

sampleNames(raw.tissue.dropped) <- c(condition.names.dropped)


##### Quality Assessment with Dropped Arrays #####

plot.colors.dropped <- rainbow(as.integer(dim(raw.tissue.dropped)[2]))

dir.create(paste(subdir.drop, "QA Images of Raw Data", sep="/"))
subdir.dropQA <- paste(subdir.drop, "QA Images of Raw Data/", sep="/")


### Boxplot ###

boxplot(raw.tissue.dropped, range=1.5, col=plot.colors.dropped, xlab="Array", ylab="Log Probe Intensity", main=paste(TISSUE, "Raw Log Probe Intensity with Dropped Arrays", sep=" "))
quartz.save(paste(subdir.dropQA, paste(TISSUE, "Raw Boxplot with Dropped Arrays.pdf", sep=" "), sep=""), type="pdf", width=15, height=7)


### Density Plot ###

hist(raw.tissue.dropped, col=plot.colors.dropped, lty=1, xlab="Log Intensity", ylab="Density", main=paste(TISSUE, "Raw Density Estimation with Dropped Arrays", sep=" "))
legend("topright", inset=0.01, cex=0.75, c(condition.names.dropped), col=plot.colors.dropped, lty=1)
quartz.save(paste(subdir.dropQA, paste(TISSUE, "Raw Density Estimation Plot with Dropped Arrays.pdf", sep=" "), sep=""), type="pdf", width=10, height=7)


### Principal Component Analysis Plot ###

pca.colors.dropped <- c(rep("red", length(LC.remaining)), rep("green", length(ROS.remaining)), rep("cyan", length(RSV.remaining)), rep("purple", length(ZC.remaining)))
pca.numbers.dropped <- c(LC.remaining, ROS.remaining, RSV.remaining, ZC.remaining)

raw.expression.matrix.dropped <- exprs(raw.tissue.dropped)
transposed.raw.expression.matrix.dropped <- t(raw.expression.matrix.dropped)

pca.values.raw.dropped <- prcomp(transposed.raw.expression.matrix.dropped)

par(mar=c(5.1, 4.1, 4.1, 6.1), xpd=T)
plot(pca.values.raw.dropped$x, col=pca.colors.dropped, pch=20, main=paste(TISSUE, "Raw PCA Plot with Dropped Arrays", sep=" "))
text(pca.values.raw.dropped$x, pos=3, offset=0.2, labels=pca.numbers.dropped, cex=0.5)
legend("topright", inset=c(-0.15,0), c(pca.conditions), cex=0.75, col=pca.legend.colors, pch=20)
quartz.save(paste(subdir.dropQA, paste(TISSUE, "Raw PCA Plot with Dropped Arrays.pdf", sep=" "), sep=""), type="pdf")
par(normal)

pca.summary.raw.dropped <- summary(pca.values.raw.dropped)
proportion.variance.raw.dropped <- pca.summary.raw.dropped$importance[2:3,1:5]

barplot(proportion.variance.raw.dropped, beside=T, col=c("black","gray"), main=paste(TISSUE, "Raw Proportion of Variance of Principal Components \n with Dropped Arrays", sep=" "), xlab="Principal Components", ylab="Percentage")
legend("topleft", inset=0.01, cex=0.75, c("Proportion of Variance", "Cumulative Proportion"), pch=15, col=c("black","gray"))
box()
quartz.save(paste(subdir.dropQA, paste(TISSUE, "Raw Proportion of Variance PCA with Dropped Arrays.pdf", sep=" "), sep=""), type="pdf", width=7, height=7)


### Hierarchical Clustering Dendogram ###

transposed.dropped <- t(raw.expression.matrix.dropped)
distance.dropped <- dist(transposed.dropped)
sample.clusters.dropped <- hclust(distance.dropped)

plot(sample.clusters.dropped, main=paste(TISSUE, "Raw Hierarchical Cluster Dendogram \n with Dropped Arrays", sep=" "), xlab="Samples", sub="")
quartz.save(paste(subdir.dropQA, paste(TISSUE, "Raw Hierarchical Cluster Dendogram with Dropped Arrays.pdf", sep=" "), sep=""), type="pdf")


### PLM Fit ###

raw.tissue.dropped.plm <- fitProbeLevelModel(raw.tissue.dropped)

dir.create(paste(subdir.dropQA, "PLM Plots", sep="/"))
subdir.dropQA.plm <- paste(subdir.dropQA, "PLM Plots/", sep="/")


# NUSE Plots #

NUSE(raw.tissue.dropped.plm, xlab="Array", main=paste(TISSUE, "Normalized Unscaled Standard Errors with Dropped Arrays", sep=" "))
quartz.save(paste(subdir.dropQA.plm, paste(TISSUE, "NUSE plot with Dropped Arrays.pdf", sep=" "), sep=""), type="pdf", width=15, height=7)


# RLE Plots #

RLE(raw.tissue.dropped.plm, xlab="Array", main=paste(TISSUE, "Relative Log Expression with Dropped Arrays", sep=" "))
quartz.save(paste(subdir.dropQA.plm, paste(TISSUE, "RLE plot with Dropped Arrays.pdf", sep=" "), sep=""), type="pdf", width=15, height=7)


##### Preprocess with Dropped Arrays #####

normalized.tissue.dropped <- rma(raw.tissue.dropped)

dir.create(paste(subdir.drop, "Preprocessed Data", sep="/"))
subdir.drop.preproc <- paste(subdir.drop, "Preprocessed Data/", sep="/")
dir.create(paste(subdir.drop.preproc, "Images", sep="/"))
subdir.drop.preproc.im <- paste(subdir.drop.preproc, "Images/", sep="/")


### Boxplot ###

boxplot(normalized.tissue.dropped, range=1.5, col=plot.colors.dropped, xlab="Array", ylab="Log Probe Intensity", main=paste(TISSUE, "Preprocessed Log Probe Intensity with Dropped Arrays", sep=" "))
quartz.save(paste(subdir.drop.preproc.im, paste(TISSUE, "Preprossed Boxplot with Dropped Arrays.pdf", sep=" "), sep=""), type="pdf", width=15, height=7)


### Density Plot ###

hist(normalized.tissue.dropped, col=plot.colors.dropped, lty=1, xlab="Log Intensity", ylab="Density", main=paste(TISSUE, "Preprocessed Density Estimation with Dropped Arrays", sep=" "))
legend("topright", inset=0.01, cex=0.75, c(condition.names.dropped), col=plot.colors.dropped, lty=1)
quartz.save(paste(subdir.drop.preproc.im, paste(TISSUE, "Preprocessed Density Estimation Plot with Dropped Arrays.pdf", sep=" "), sep=""), type="pdf", width=10, height=7)


### Principal Component Analysis ###

preprocessed.expression.matrix.dropped <- exprs(normalized.tissue.dropped)
transposed.preprocessed.expression.matrix <- t(preprocessed.expression.matrix.dropped)

pca.values.preprocessed.dropped <- prcomp(transposed.preprocessed.expression.matrix)

par(mar=c(5.1, 4.1, 4.1, 6.1), xpd=T)
plot(pca.values.preprocessed.dropped$x, col=pca.colors.dropped, pch=20, main=paste(TISSUE, "Preprocessed PCA Plot with Dropped Arrays", sep=" "))
text(pca.values.preprocessed.dropped$x, pos=3, offset=0.2, labels=pca.numbers.dropped, cex=0.5)
legend("topright", inset=c(-0.15,0), c(pca.conditions), cex=0.75, col=pca.legend.colors, pch=20)
quartz.save(paste(subdir.drop.preproc.im, paste(TISSUE, "Preprocessed PCA Plot with Dropped Arrays.pdf", sep=" "), sep=""), type="pdf")
par(normal)

pca.summary.preprocessed.dropped <- summary(pca.values.preprocessed.dropped)
proportion.variance.preprocessed.dropped <- pca.summary.preprocessed.dropped$importance[2:3,1:5]

barplot(proportion.variance.preprocessed.dropped, beside=T, col=c("black","gray"), main=paste(TISSUE, "Preprocessed Proportion of Variance of Principal Components \n with Dropped Arrays", sep=" "), xlab="Principal Components", ylab="Percentage")
legend("topleft", inset=0.01, cex=0.75, c("Proportion of Variance", "Cumulative Proportion"), pch=15, col=c("black","gray"))
box()
quartz.save(paste(subdir.drop.preproc.im, paste(TISSUE, "Preprocessed Proportion of Variance PCA with Dropped Arrays.pdf", sep=" "), sep=""), type="pdf", width=7, height=7)


### Hierarchical Clustering Dendogram ###

preprocessed.transposed.dropped <- t(preprocessed.expression.matrix.dropped)
preprocessed.distance.dropped <- dist(preprocessed.transposed.dropped)
preprocessed.sample.clusters.dropped <- hclust(preprocessed.distance.dropped)

plot(preprocessed.sample.clusters.dropped, main=paste(TISSUE, "Preprocessed Hierarchical Cluster Dendogram \n with Dropped Arrays", sep=" "), xlab="Samples", sub="")
quartz.save(paste(subdir.drop.preproc.im, paste(TISSUE, "Preprocessed Hierarchical Cluster Dendogram with Dropped Arrays.pdf", sep=" "), sep=""), type="pdf")


##### Filter Genes to Remove Those with No Annotation #####

dir.create(paste(subdir.drop.preproc, "Filtered Genes Data", sep="/"))
subdir.drop.preproc.filt <- paste(subdir.drop.preproc, "Filtered Genes Data/", sep="/")

annotation.file <- read.csv("~/Desktop/Annotations for Rat Gene 21st.csv", header=TRUE)
annotation.matrix <- as.matrix(annotation.file)
annotations <- annotation.matrix[order(annotation.file[,1], annotation.file[,2]),]
#creates a matrix where the genes are ordered by their transcript cluster

annotated.values <- cbind(annotations, preprocessed.expression.matrix.dropped)
#combines gene annotations with expression values, which are already ordered according to transcript cluster
annotated.gene.values.all <- annotated.values[,-1]
#removes transcript cluster column
annotated.gene.values.no.NA <- na.omit(annotated.gene.values.all)
#removes all genes with no annotation

rownames(annotated.gene.values.no.NA) <- annotated.gene.values.no.NA[,1]
annotated.gene.values <- annotated.gene.values.no.NA[,-1]
#removes column of gene IDs so all columns are expression values


##### Average Expression Values for Duplicated Genes #####

genes <- as.numeric(rownames(annotated.gene.values))
#warning is ok

duplicates <- data.frame(Gene=genes, annotated.gene.values)
#warning is ok

remove.duplicates <- aggregate(duplicates, by=list(duplicates[,1]), FUN=mean)
#for genes that are duplicated, for each column, the values are averaged

rownames(remove.duplicates)<-remove.duplicates[,2]

no.duplicates <- remove.duplicates[,-c(1:2)]


### Create Reference List of All Annotated Genes for FunNet ##

reference.list <- as.numeric(rownames(no.duplicates))

write.table(reference.list, paste(subdir.drop.preproc.filt, paste(TISSUE, "Reference Gene List.txt", sep=" "), sep=""), quote=F, row.names=F, col.names=F)


##### Images for Filtered Genes Data #####

dir.create(paste(subdir.drop.preproc.filt, "Images", sep="/"))
subsubdir.images <- paste(subdir.drop.preproc.filt, "Images/", sep="/")


### Boxplot ###

boxplot(no.duplicates, range=1.5, col=plot.colors.dropped, xlab="Array", ylab="Log Probe Intensity", main=paste(TISSUE, "Preprocessed Log Probe Intensity with Annotated Genes", sep=" "))
quartz.save(paste(subsubdir.images, paste(TISSUE, "Preprossed Boxplot with Annotated Genes.pdf", sep=" "), sep=""), type="pdf", width=15, height=7)


### PCA Plot of Filtered Genes ###

transposed.annotated.gene.values <- t(no.duplicates)
pca.annotated.genes <- prcomp(transposed.annotated.gene.values)

par(mar=c(5.1, 4.1, 4.1, 6.1), xpd=T)
plot(pca.annotated.genes$x, col=pca.colors.dropped, pch=20, main=paste(TISSUE, "Preprocessed PCA Plot with Annotated Genes", sep=" "))
text(pca.annotated.genes$x, pos=3, offset=0.2, labels=pca.numbers.dropped, cex=0.5)
legend("topright", inset=c(-0.15,0), c(pca.conditions), cex=0.75, col=pca.legend.colors, pch=20)
quartz.save(paste(subsubdir.images, paste(TISSUE, "Preprocessed PCA Plot with Annotated Genes", sep=" "), sep=""), type="pdf")
par(normal)

pca.summary <- summary(pca.annotated.genes)
proportion.variance <- pca.summary$importance[2:3,1:5]

barplot(proportion.variance,beside=T, col=c("black","gray"), main=paste(TISSUE, "Preprocessed Proportion of Variance of Principal Components \n with Annotated Genes", sep=" "), xlab="Principal Components", ylab="Percentage")
legend("topleft",inset=0.01, cex=0.75, c("Proportion of Variance", "Cumulative Proportion"), pch=15, col=c("black","gray"))
box()
quartz.save(paste(subsubdir.images, paste(TISSUE, "Preprocessed Proportion of Variance PCA with Annotated Genes.pdf", sep=" "), sep=""), type="pdf", width=7, height=7)


### Hierarchical Clustering Dendogram ###

annotated.transposed <- t(no.duplicates)
annotated.distance <- dist(annotated.transposed)
annotated.sample.clusters <- hclust(annotated.distance)

plot(annotated.sample.clusters, main=paste(TISSUE, "Preprocessed Hierarchical Cluster Dendogram \n with Annotated Genes", sep=" "), xlab="Samples", sub="")
quartz.save(paste(subsubdir.images, paste(TISSUE, "Preprocessed Hierarchical Cluster Dendogram with Annotated Genes.pdf", sep=" "), sep=""), type="pdf")


##### Use ANOVA to Determine Differentially Expressed Genes #####

#code modified from Pavlidis, P. 2003 Methods 31(4):282-289

dir.create(paste(subdir.drop.preproc.filt, "Differentially Expressed Gene Lists", sep="/"))
subsubdir.DE <- paste(subdir.drop.preproc.filt, "Differentially Expressed Gene Lists/", sep="/")

remaining.conditions <- c(rep("LC", length(LC.remaining)), rep("ROS", length(ROS.remaining)), rep("RSV", length(RSV.remaining)), rep("ZC", length(ZC.remaining)))

conditions <- factor(remaining.conditions)
#allows you to perform one-way ANOVA with unequal sample sizes

anova.function <- function(x){
	data <- data.frame(conditions,x)
	anova(aov(x~conditions, data))
}

anova.results <- apply(no.duplicates, 1, anova.function)
#applies ANOVA function to each row (gene) of the filtered genes

pvalue.function <- function(x){
	x["Pr(>F)"][1,]
}

pvalues <- data.frame(lapply(anova.results, pvalue.function))

pvalues.table <- t(pvalues)
colnames(pvalues.table) <- "Condition P Value"


### False-Discovery Rate (FDR) Correction ###

adjusted.pvalues <- as.data.frame(p.adjust(pvalues.table[,1], method="fdr"))
colnames(adjusted.pvalues) <- "Adjusted P Value"

p.05 <- which(adjusted.pvalues[,1]<0.05)

diff.expressed.genes <- no.duplicates[p.05,]
#gets a list of genes that have a p-value of less than 0.05 -- differentially expressed

ANOVA.geneIDs <- rownames(diff.expressed.genes)

write.table(ANOVA.geneIDs, paste(subsubdir.DE, paste(TISSUE, "ANOVA DE Genes.txt", sep=" "), sep=""), quote=F, row.names=F, col.names=F)


##### Post-Hoc Pairwise Comparisons with FDR #####

dir.create(paste(subsubdir.DE, "Pairwise DE Genes", sep="/"))
subsubdir.DE.pair <- paste(subsubdir.DE, "Pairwise DE Genes/", sep="/")

pairwise.comparisons <- c("ROSvLC", "RSVvLC", "RSVvROS", "ZCvLC", "ROSvZC", "RSVvZC")

diff.expressed.matrix <- as.matrix(diff.expressed.genes)

ttest.results <- matrix(nrow=nrow(diff.expressed.matrix), ncol=6)
colnames(ttest.results) <- pairwise.comparisons

for (i in 1:nrow(diff.expressed.matrix)){
	p <- pairwise.t.test(diff.expressed.matrix[i,], conditions, p.adj="fdr")
	pv <- p$p.value
	ttest.results[i,1] <- pv[1,1]
	ttest.results[i,2] <- pv[2,1]
	ttest.results[i,3] <- pv[2,2]
	ttest.results[i,4] <- pv[3,1]
	ttest.results[i,5] <- pv[3,2]
	ttest.results[i,6] <- pv[3,3]
}

rownames(ttest.results) <- rownames(diff.expressed.genes)

pairwise.diff.genes <- matrix(nrow=nrow(diff.expressed.genes), ncol=6)

for (i in 1:nrow(ttest.results)){
	if (ttest.results[i,1]<=0.05){
		pairwise.diff.genes[i,1] <- (i)
	}
	if (ttest.results[i,2]<=0.05){
		pairwise.diff.genes[i,2] <- (i)
	}
	if (ttest.results[i,3]<=0.05){
		pairwise.diff.genes[i,3] <- (i)
	}
	if (ttest.results[i,4]<=0.05){
		pairwise.diff.genes[i,4] <- (i)
	}
	if (ttest.results[i,5]<=0.05){
		pairwise.diff.genes[i,5] <- (i)
	}
	if (ttest.results[i,6]<=0.05){
		pairwise.diff.genes[i,6] <- (i)
	}
}

ROSvLC.info <- na.omit(pairwise.diff.genes[,1])
RSVvLC.info <- na.omit(pairwise.diff.genes[,2])
RSVvROS.info <- na.omit(pairwise.diff.genes[,3])
ZCvLC.info <- na.omit(pairwise.diff.genes[,4])
ROSvZC.info <- na.omit(pairwise.diff.genes[,5])
RSVvZC.info <- na.omit(pairwise.diff.genes[,6])
#gets the row numbers for the differentially expressed genes for each pairwise comparison

ROSvLC.diff.values <- as.matrix(ttest.results[,1][ROSvLC.info])
RSVvLC.diff.values <- as.matrix(ttest.results[,2][RSVvLC.info])
RSVvROS.diff.values <- as.matrix(ttest.results[,3][RSVvROS.info])
ZCvLC.diff.values <- as.matrix(ttest.results[,4][ZCvLC.info])
ROSvZC.diff.values <- as.matrix(ttest.results[,5][ROSvZC.info])
RSVvZC.diff.values <- as.matrix(ttest.results[,6][RSVvZC.info])
#gets the expression values for the differentially expressed genes for each pairwise comparison

ROSvLC.genes <- as.numeric(rownames(ROSvLC.diff.values))
RSVvLC.genes <- as.numeric(rownames(RSVvLC.diff.values))
RSVvROS.genes <- as.numeric(rownames(RSVvROS.diff.values))
ZCvLC.genes <- as.numeric(rownames(ZCvLC.diff.values))
ROSvZC.genes <- as.numeric(rownames(ROSvZC.diff.values))
RSVvZC.genes <- as.numeric(rownames(RSVvZC.diff.values))
#gets a list of just the differentially expressed gene IDs for each pairwise comparison


### Save Gene Lists to Files ###

dir.create(paste(subsubdir.DE.pair, "All DE Genes", sep="/"))
subsubdir.DE.pairA<-paste(subsubdir.DE.pair,"All DE Genes/",sep="/")

write.table(ROSvLC.genes, paste(subsubdir.DE.pairA, paste(TISSUE, "ROSvLC DE Genes.txt", sep=" "), sep=""), quote=F, row.names=F, col.names=F)
write.table(RSVvLC.genes, paste(subsubdir.DE.pairA, paste(TISSUE, "RSVvLC DE Genes.txt", sep=" "), sep=""), quote=F, row.names=F, col.names=F)
write.table(RSVvROS.genes,paste(subsubdir.DE.pairA, paste(TISSUE, "RSVvROS DE Genes.txt", sep=" "), sep=""), quote=F, row.names=F, col.names=F)
write.table(ZCvLC.genes,paste(subsubdir.DE.pairA, paste(TISSUE, "ZCvLC DE Genes.txt", sep=" "), sep=""), quote=F, row.names=F, col.names=F)
write.table(ROSvZC.genes,paste(subsubdir.DE.pairA, paste(TISSUE, "ROSvZC DE Genes.txt", sep=" "), sep=""), quote=F, row.names=F, col.names=F)
write.table(RSVvZC.genes, paste(subsubdir.DE.pairA, paste(TISSUE, "RSVvZC DE Genes.txt", sep=" "), sep=""), quote=F, row.names=F, col.names=F)


##### Determine Up- and Down-Regulated Genes #####

average.condition.names <- c(rep("LC", length(LC.remaining)), rep("ROS", length(ROS.remaining)), rep("RSV", length(RSV.remaining)), rep("ZC", length(ZC.remaining)))

condition.value.columns <- diff.expressed.matrix

colnames(condition.value.columns) <- average.condition.names
#replaces the column names representing each individual array number with general names representing the conditions, so that the expression values for each condition can be averaged (next steps)

transposed.condition.values <- t(condition.value.columns)

group.condition.values <- data.frame(Condition= average.condition.names,  transposed.condition.values)
#warning is ok

mean.condition.values <- aggregate(group.condition.values, by=list(group.condition.values[,1]), FUN=mean)
#warning is ok
#averages out the values for each condition for each gene

combined.condition.values <- t(mean.condition.values)
cvalues <- data.frame(combined.condition.values[-c(1,2),])

LC.values <- as.numeric(levels(cvalues$X1))[cvalues$X1]
ROS.values <- as.numeric(levels(cvalues$X2))[cvalues$X2]
RSV.values <- as.numeric(levels(cvalues$X3))[cvalues$X3]
ZC.values <- as.numeric(levels(cvalues$X4))[cvalues$X4]
#the aggragate function results in the columns being identified as factors instead of numbers; this changes each number from a level of a factor to a numeric value

average.condition.values <- data.frame(LC=LC.values, ROS=ROS.values, RSV=RSV.values, ZC=ZC.values)

gene.names <- rownames(diff.expressed.genes)

rownames(average.condition.values) <- gene.names

pairwise.differences <- data.frame(Gene=as.numeric(gene.names), ROSvLC=(average.condition.values[,2]-average.condition.values[,1]), RSVvLC=(average.condition.values[,3]-average.condition.values[,1]), RSVvROS=(average.condition.values[,3]-average.condition.values[,2]), ZCvLC=(average.condition.values[,4]-average.condition.values[,1]), ROSvZC=(average.condition.values[,2]-average.condition.values[,4]), RSVvZC=(average.condition.values[,3]-average.condition.values[,4]))
#subtracts the expression values of the second condition from the first condition for each pairing, to determine whether each gene was up- or down-regulated

ROSvLC.difference <- pairwise.differences[,1:2]
RSVvLC.difference <- pairwise.differences[,c(1,3)]
RSVvROS.difference <- pairwise.differences[,c(1,4)]
ZCvLC.difference <- pairwise.differences[,c(1,5)]
ROSvZC.difference <- pairwise.differences[,c(1,6)]
RSVvZC.difference <- pairwise.differences[,c(1,7)]
#combines the gene ID with the difference (positive or negative) for all ANOVA differentially expressed genes for each pairwise comparison

list.ROSvLC <- ROSvLC.difference[ROSvLC.info,]
list.RSVvLC <- RSVvLC.difference[RSVvLC.info,]
list.RSVvROS <- RSVvROS.difference[RSVvROS.info,]
list.ZCvLC <- ZCvLC.difference[ZCvLC.info,]
list.ROSvZC <- ROSvZC.difference[ROSvZC.info,]
list.RSVvZC <- RSVvZC.difference[RSVvZC.info,]
#reduces the genes to just those determined to be differentially expressed in the pairwise conditions by t tests

ROSvLC.updown <- data.frame(Upregulated=(NA) ,Downregulated=(NA))
RSVvLC.updown <- data.frame(Upregulated=(NA), Downregulated=(NA))
RSVvROS.updown <- data.frame(Upregulated=(NA), Downregulated=(NA))
ZCvLC.updown <- data.frame(Upregulated=(NA), Downregulated=(NA))
ROSvZC.updown <- data.frame(Upregulated=(NA), Downregulated=(NA))
RSVvZC.updown <- data.frame(Upregulated=(NA), Downregulated=(NA))

for (i in 1:nrow(list.ROSvLC)){
	if (list.ROSvLC[i,2]>0){
		ROSvLC.updown[i,1] <- list.ROSvLC[i,1]
	}
	if (list.ROSvLC[i,2]<0){
		ROSvLC.updown[i,2] <- list.ROSvLC[i,1]
	}
}

for (i in 1:nrow(list.RSVvLC)){
	if (list.RSVvLC[i,2]>0){
		RSVvLC.updown[i,1] <- list.RSVvLC[i,1]
	}
	if (list.RSVvLC[i,2]<0){
		RSVvLC.updown[i,2] <- list.RSVvLC[i,1]
	}
}

for (i in 1:nrow(list.RSVvROS)){
	if (list.RSVvROS[i,2]>0){
		RSVvROS.updown[i,1] <- list.RSVvROS[i,1]
	}
	if (list.RSVvROS[i,2]<0){
		RSVvROS.updown[i,2] <- list.RSVvROS[i,1]
	}
}

for (i in 1:nrow(list.ZCvLC)){
	if (list.ZCvLC[i,2]>0){
		ZCvLC.updown[i,1] <- list.ZCvLC[i,1]
	}
	if (list.ZCvLC[i,2]<0){
		ZCvLC.updown[i,2] <- list.ZCvLC[i,1]
	}
}

for (i in 1:nrow(list.ROSvZC)){
	if (list.ROSvZC[i,2]>0){
		ROSvZC.updown[i,1] <- list.ROSvZC[i,1]
	}
	if (list.ROSvZC[i,2]<0){
		ROSvZC.updown[i,2] <- list.ROSvZC[i,1]
	}
}

for (i in 1:nrow(list.RSVvZC)){
	if (list.RSVvZC[i,2]>0){
		RSVvZC.updown[i,1] <- list.RSVvZC[i,1]
	}
	if (list.RSVvZC[i,2]<0){
		RSVvZC.updown[i,2] <- list.RSVvZC[i,1]
	}
}
#puts gene ID into the up- or down-regulated column appropriately for each pariwise comparison

ROSvLC.up <- na.omit(ROSvLC.updown[,1])
RSVvLC.up <- na.omit(RSVvLC.updown[,1])
RSVvROS.up <- na.omit(RSVvROS.updown[,1])
ZCvLC.up <- na.omit(ZCvLC.updown[,1])
ROSvZC.up <- na.omit(ROSvZC.updown[,1])
RSVvZC.up <- na.omit(RSVvZC.updown[,1])

ROSvLC.down <- na.omit(ROSvLC.updown[,2])
RSVvLC.down <- na.omit(RSVvLC.updown[,2])
RSVvROS.down <- na.omit(RSVvROS.updown[,2])
ZCvLC.down <- na.omit(ZCvLC.updown[,2])
ROSvZC.down <- na.omit(ROSvZC.updown[,2])
RSVvZC.down <- na.omit(RSVvZC.updown[,2])
#gets rid of NA values and separates the up- and down-regulated gene IDs into separate lists


### Write Up- and Down-Regulated Genes to Files for FunNet ###

dir.create(paste(subsubdir.DE.pair, "Up Regulated Genes", sep="/"))
subsubdir.DE.pairU <- paste(subsubdir.DE.pair, "Up Regulated Genes/", sep="/")

write.table(ROSvLC.up, paste(subsubdir.DE.pairU, paste(TISSUE, "ROSvLC Upregulated Genes.txt", sep=" "), sep=""), quote=F, row.names=F, col.names=F)
write.table(RSVvLC.up, paste(subsubdir.DE.pairU, paste(TISSUE, "RSVvLC Upregulated Genes.txt", sep=" "), sep=""), quote=F, row.names=F, col.names=F)
write.table(RSVvROS.up, paste(subsubdir.DE.pairU, paste(TISSUE, "RSVvROS Upregulated Genes.txt", sep=" "), sep=""), quote=F, row.names=F, col.names=F)
write.table(ZCvLC.up, paste(subsubdir.DE.pairU, paste(TISSUE, "ZCvLC Upregulated Genes.txt", sep=" "), sep=""), quote=F, row.names=F, col.names=F)
write.table(ROSvZC.up, paste(subsubdir.DE.pairU, paste(TISSUE, "ROSvZC Upregulated Genes.txt", sep=" "), sep=""), quote=F, row.names=F, col.names=F)
write.table(RSVvZC.up, paste(subsubdir.DE.pairU, paste(TISSUE, "RSVvZC Upregulated Genes.txt", sep=" "), sep=""), quote=F, row.names=F, col.names=F)

dir.create(paste(subsubdir.DE.pair, "Down Regulated Genes", sep="/"))
subsubdir.DE.pairD <- paste(subsubdir.DE.pair, "Down Regulated Genes/", sep="/")

write.table(ROSvLC.down, paste(subsubdir.DE.pairD, paste(TISSUE, "ROSvLC Downregulated Genes.txt", sep=" "), sep=""), quote=F, row.names=F, col.names=F)
write.table(RSVvLC.down, paste(subsubdir.DE.pairD, paste(TISSUE, "RSVvLC Downregulated Genes.txt", sep=" "), sep=""), quote=F, row.names=F, col.names=F)
write.table(RSVvROS.down, paste(subsubdir.DE.pairD, paste(TISSUE, "RSVvROS Downregulated Genes.txt", sep=" "), sep=""), quote=F, row.names=F, col.names=F)
write.table(ZCvLC.down, paste(subsubdir.DE.pairD, paste(TISSUE, "ZCvLC Downregulated Genes.txt", sep=" "), sep=""), quote=F, row.names=F, col.names=F)
write.table(ROSvZC.down, paste(subsubdir.DE.pairD, paste(TISSUE, "ROSvZC Downregulated Genes.txt", sep=" "), sep=""), quote=F, row.names=F, col.names=F)
write.table(RSVvZC.down, paste(subsubdir.DE.pairD, paste(TISSUE, "RSVvZC Downregulated Genes.txt", sep=" "), sep=""), quote=F, row.names=F, col.names=F)


##### Determine Number of DE, Up- and Down-Regulated Genes #####

total.num <- c(length(ROSvLC.genes), length(RSVvLC.genes), length(RSVvROS.genes), length(ZCvLC.genes), length(ROSvZC.genes), length(RSVvZC.genes))
up.num <- c(length(ROSvLC.up), length(RSVvLC.up), length(RSVvROS.up), length(ZCvLC.up), length(ROSvZC.up), length(RSVvZC.up))
down.num <- c(length(ROSvLC.down), length(RSVvLC.down), length(RSVvROS.down), length(ZCvLC.down), length(ROSvZC.down), length(RSVvZC.down))

DE.numbers <- data.frame(Comparison=(pairwise.comparisons), Total=(total.num), UpRegulated=(up.num), DownRegulated=(down.num))

write.table(DE.numbers, paste(subsubdir.DE, paste(TISSUE, "Differentially Expressed Gene Numbers.txt", sep=" "), sep=""), quote=F, row.names=F, sep="\t")

ANOVA.DE.number <- c("Total DE Genes Determined by ANOVA",length(p.05))

write.table(c("\n", ANOVA.DE.number), paste(subsubdir.DE, paste(TISSUE, "Differentially Expressed Gene Numbers.txt", sep=" "), sep=""), append=T, quote=F, row.names=F, col.names=F, sep="\t")
