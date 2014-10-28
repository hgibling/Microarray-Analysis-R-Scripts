library(plyr)
library(limma)
library(venneuler)

##### Determine Overlapped Genes for Rosiglitazone #####

### Total ###

ROSvLC.id <- cbind(ROSvLC.genes, 1)
ROSvZC.id <- cbind(ROSvZC.genes, 3)
ZCvLC.id <- cbind(ZCvLC.genes, 5)

ROS.genes <- as.data.frame(rbind(ROSvLC.id, ROSvZC.id, ZCvLC.id))

colnames(ROS.genes) <- c("GeneID", "Pair")

ROS.ids <- ddply(ROS.genes, "GeneID", summarize, Pair=sum(Pair))

ROS.LC <- which(ROS.ids[,2]==1)
ROS.ZC <- which(ROS.ids[,2]==3)
ROS.ZCLC <- which(ROS.ids[,2]==5)
ROS.LC.ZC <- which(ROS.ids[,2]==4)
ROS.LC.ZCLC <- which(ROS.ids[,2]==6)
ROS.ZC.ZCLC <- which(ROS.ids[,2]==8)
ROS.all <- which(ROS.ids[,2]==9)

ROS.matrix <- matrix(nrow=nrow(ROS.ids), ncol=3)

for (i in 1:nrow(ROS.ids)){
	if (ROS.ids[i,2]==1){
		ROS.matrix[i,1] <- 1
		ROS.matrix[i,2:3] <- 0
	}
	if (ROS.ids[i,2]==3){
		ROS.matrix[i,2] <- 1
		ROS.matrix[i,c(1,3)] <- 0
	}
	if (ROS.ids[i,2]==5){
		ROS.matrix[i,3] <- 1
		ROS.matrix[i,1:2] <- 0
	}
	if (ROS.ids[i,2]==4){
		ROS.matrix[i,1:2] <- 1
		ROS.matrix[i,3] <- 0
	}
	if (ROS.ids[i,2]==6){
		ROS.matrix[i,c(1,3)] <- 1
		ROS.matrix[i,2] <- 0
	}
	if (ROS.ids[i,2]==8){
		ROS.matrix[i,2:3] <- 1
		ROS.matrix[i,1] <- 0
	}
	if (ROS.ids[i,2]==9){
		ROS.matrix[i,1:3] <- 1
	}
}


### Up-Regulated ###

ROSvLC.id.up <- cbind(ROSvLC.up, 1)
ROSvZC.id.up <- cbind(ROSvZC.up, 3)
ZCvLC.id.up <- cbind(ZCvLC.up, 5)

ROS.genes.up <- as.data.frame(rbind(ROSvLC.id.up, ROSvZC.id.up, ZCvLC.id.up))

colnames(ROS.genes.up) <- c("GeneID", "Pair")

ROS.ids.up <- ddply(ROS.genes.up, "GeneID", summarize, Pair=sum(Pair))

ROS.LC.up <- which(ROS.ids.up[,2]==1)
ROS.ZC.up <- which(ROS.ids.up[,2]==3)
ROS.ZCLC.up <- which(ROS.ids.up[,2]==5)
ROS.LC.ZC.up <- which(ROS.ids.up[,2]==4)
ROS.LC.ZCLC.up <- which(ROS.ids.up[,2]==6)
ROS.ZC.ZCLC.up <- which(ROS.ids.up[,2]==8)
ROS.all.up <- which(ROS.ids.up[,2]==9)

ROS.matrix.up <- matrix(nrow=nrow(ROS.ids.up), ncol=3)

for (i in 1:nrow(ROS.ids.up)){
	if (ROS.ids.up[i,2]==1){
		ROS.matrix.up[i,1] <- 1
		ROS.matrix.up[i,2:3] <- 0
	}
	if (ROS.ids.up[i,2]==3){
		ROS.matrix.up[i,2] <- 1
		ROS.matrix.up[i,c(1,3)] <- 0
	}
	if (ROS.ids.up[i,2]==5){
		ROS.matrix.up[i,3] <- 1
		ROS.matrix.up[i,1:2] <- 0
	}
	if (ROS.ids.up[i,2]==4){
		ROS.matrix.up[i,1:2] <- 1
		ROS.matrix.up[i,3] <- 0
	}
	if (ROS.ids.up[i,2]==6){
		ROS.matrix.up[i,c(1,3)] <- 1
		ROS.matrix.up[i,2] <- 0
	}
	if (ROS.ids.up[i,2]==8){
		ROS.matrix.up[i,2:3] <- 1
		ROS.matrix.up[i,1] <- 0
	}
	if (ROS.ids.up[i,2]==9){
		ROS.matrix.up[i,1:3] <- 1
	}
}


### Down-Regulated ###

ROSvLC.id.down <- cbind(ROSvLC.down, 1)
ROSvZC.id.down <- cbind(ROSvZC.down, 3)
ZCvLC.id.down <- cbind(ZCvLC.down, 5)

ROS.genes.down <- as.data.frame(rbind(ROSvLC.id.down, ROSvZC.id.down, ZCvLC.id.down))

colnames(ROS.genes.down) <- c("GeneID", "Pair")

ROS.ids.down <- ddply(ROS.genes.down, "GeneID", summarize, Pair=sum(Pair))

ROS.LC.down <- which(ROS.ids.down[,2]==1)
ROS.ZC.down <- which(ROS.ids.down[,2]==3)
ROS.ZCLC.down <- which(ROS.ids.down[,2]==5)
ROS.LC.ZC.down <- which(ROS.ids.down[,2]==4)
ROS.LC.ZCLC.down <- which(ROS.ids.down[,2]==6)
ROS.ZC.ZCLC.down <- which(ROS.ids.down[,2]==8)
ROS.all.down <- which(ROS.ids.down[,2]==9)

ROS.matrix.down <- matrix(nrow=nrow(ROS.ids.down), ncol=3)

for (i in 1:nrow(ROS.ids.down)){
	if (ROS.ids.down[i,2]==1){
		ROS.matrix.down[i,1] <- -1
		ROS.matrix.down[i,2:3] <- 0
	}
	if (ROS.ids.down[i,2]==3){
		ROS.matrix.down[i,2] <- -1
		ROS.matrix.down[i,c(1,3)] <- 0
	}
	if (ROS.ids.down[i,2]==5){
		ROS.matrix.down[i,3] <- -1
		ROS.matrix.down[i,1:2] <- 0
	}
	if (ROS.ids.down[i,2]==4){
		ROS.matrix.down[i,1:2] <- -1
		ROS.matrix.down[i,3] <- 0
	}
	if (ROS.ids.down[i,2]==6){
		ROS.matrix.down[i,c(1,3)] <- -1
		ROS.matrix.down[i,2] <- 0
	}
	if (ROS.ids.down[i,2]==8){
		ROS.matrix.down[i,2:3] <- -1
		ROS.matrix.down[i,1] <- 0
	}
	if (ROS.ids.down[i,2]==9){
		ROS.matrix.down[i,1:3] <- -1
	}
}

ROS.matrix.updown <- rbind(ROS.matrix.up, ROS.matrix.down)


### Venn Diagrams ###

par(mfrow=c(1,2))
vennDiagram(ROS.matrix, names=c("ROSvLC", "ROSvZC", "ZCvLC"), circle.col=c("darkgreen", "orange", "blue"), main=paste(TISSUE, "Overlap in All DE Genes", sep=" "))
vennDiagram(ROS.matrix.updown, include=c("up", "down"), names=c("ROSvLC", "ROSvZC", "ZCvLC"), circle.col=c("darkgreen", "orange", "blue"), main=paste(TISSUE, "Overlap in DE Genes", sep=" "), sub="test")
quartz.save(paste("~/Desktop/", paste(TISSUE, "Venn ROS.pdf", sep=" "), sep=""), type="pdf", width=13, height=7)


##### Determine Overlapped Genes for Resveratrol #####

### Total ###

RSVvLC.id <- cbind(RSVvLC.genes, 1)
RSVvZC.id <- cbind(RSVvZC.genes, 3)

RSV.genes <- as.data.frame(rbind(RSVvLC.id, RSVvZC.id, ZCvLC.id))

colnames(RSV.genes) <- c("GeneID", "Pair")

RSV.ids <- ddply(RSV.genes, "GeneID", summarize, Pair=sum(Pair))

RSV.LC <- which(RSV.ids[,2]==1)
RSV.ZC <- which(RSV.ids[,2]==3)
RSV.ZCLC <- which(RSV.ids[,2]==5)
RSV.LC.ZC <- which(RSV.ids[,2]==4)
RSV.LC.ZCLC <- which(RSV.ids[,2]==6)
RSV.ZC.ZCLC <- which(RSV.ids[,2]==8)
RSV.all <- which(RSV.ids[,2]==9)

RSV.matrix <- matrix(nrow=nrow(RSV.ids), ncol=3)

for (i in 1:nrow(RSV.ids)){
	if (RSV.ids[i,2]==1){
		RSV.matrix[i,1] <- 1
		RSV.matrix[i,2:3] <- 0
	}
	if (RSV.ids[i,2]==3){
		RSV.matrix[i,2] <- 1
		RSV.matrix[i,c(1,3)] <- 0
	}
	if (RSV.ids[i,2]==5){
		RSV.matrix[i,3] <- 1
		RSV.matrix[i,1:2] <- 0
	}
	if (RSV.ids[i,2]==4){
		RSV.matrix[i,1:2] <- 1
		RSV.matrix[i,3] <- 0
	}
	if (RSV.ids[i,2]==6){
		RSV.matrix[i,c(1,3)] <- 1
		RSV.matrix[i,2] <- 0
	}
	if (RSV.ids[i,2]==8){
		RSV.matrix[i,2:3] <- 1
		RSV.matrix[i,1] <- 0
	}
	if (RSV.ids[i,2]==9){
		RSV.matrix[i,1:3] <- 1
	}
}


### Up-Regulated ###

RSVvLC.id.up <- cbind(RSVvLC.up, 1)
RSVvZC.id.up <- cbind(RSVvZC.up, 3)

RSV.genes.up <- as.data.frame(rbind(RSVvLC.id.up, RSVvZC.id.up, ZCvLC.id.up))

colnames(RSV.genes.up) <- c("GeneID", "Pair")

RSV.ids.up <- ddply(RSV.genes.up, "GeneID", summarize, Pair=sum(Pair))

RSV.LC.up <- which(RSV.ids.up[,2]==1)
RSV.ZC.up <- which(RSV.ids.up[,2]==3)
RSV.ZCLC.up <- which(RSV.ids.up[,2]==5)
RSV.LC.ZC.up <- which(RSV.ids.up[,2]==4)
RSV.LC.ZCLC.up <- which(RSV.ids.up[,2]==6)
RSV.ZC.ZCLC.up <- which(RSV.ids.up[,2]==8)
RSV.all.up <- which(RSV.ids.up[,2]==9)

RSV.matrix.up <- matrix(nrow=nrow(RSV.ids.up), ncol=3)

for (i in 1:nrow(RSV.ids.up)){
	if (RSV.ids.up[i,2]==1){
		RSV.matrix.up[i,1] <- 1
		RSV.matrix.up[i,2:3] <- 0
	}
	if (RSV.ids.up[i,2]==3){
		RSV.matrix.up[i,2] <- 1
		RSV.matrix.up[i,c(1,3)] <- 0
	}
	if (RSV.ids.up[i,2]==5){
		RSV.matrix.up[i,3] <- 1
		RSV.matrix.up[i,1:2] <- 0
	}
	if (RSV.ids.up[i,2]==4){
		RSV.matrix.up[i,1:2] <- 1
		RSV.matrix.up[i,3] <- 0
	}
	if (RSV.ids.up[i,2]==6){
		RSV.matrix.up[i,c(1,3)] <- 1
		RSV.matrix.up[i,2] <- 0
	}
	if (RSV.ids.up[i,2]==8){
		RSV.matrix.up[i,2:3] <- 1
		RSV.matrix.up[i,1] <- 0
	}
	if (RSV.ids.up[i,2]==9){
		RSV.matrix.up[i,1:3] <- 1
	}
}


### Down-Regulated ###

RSVvLC.id.down <- cbind(RSVvLC.down, 1)
RSVvZC.id.down <- cbind(RSVvZC.down, 3)

RSV.genes.down <- as.data.frame(rbind(RSVvLC.id.down, RSVvZC.id.down, ZCvLC.id.down))

colnames(RSV.genes.down) <- c("GeneID", "Pair")

RSV.ids.down <- ddply(RSV.genes.down, "GeneID", summarize, Pair=sum(Pair))

RSV.LC.down <- which(RSV.ids.down[,2]==1)
RSV.ZC.down <- which(RSV.ids.down[,2]==3)
RSV.ZCLC.down <- which(RSV.ids.down[,2]==5)
RSV.LC.ZC.down <- which(RSV.ids.down[,2]==4)
RSV.LC.ZCLC.down <- which(RSV.ids.down[,2]==6)
RSV.ZC.ZCLC.down <- which(RSV.ids.down[,2]==8)
RSV.all.down <- which(RSV.ids.down[,2]==9)

RSV.matrix.down <- matrix(nrow=nrow(RSV.ids.down), ncol=3)

for (i in 1:nrow(RSV.ids.down)){
	if (RSV.ids.down[i,2]==1){
		RSV.matrix.down[i,1] <- -1
		RSV.matrix.down[i,2:3] <- 0
	}
	if (RSV.ids.down[i,2]==3){
		RSV.matrix.down[i,2] <- -1
		RSV.matrix.down[i,c(1,3)] <- 0
	}
	if (RSV.ids.down[i,2]==5){
		RSV.matrix.down[i,3] <- -1
		RSV.matrix.down[i,1:2] <- 0
	}
	if (RSV.ids.down[i,2]==4){
		RSV.matrix.down[i,1:2] <- -1
		RSV.matrix.down[i,3] <- 0
	}
	if (RSV.ids.down[i,2]==6){
		RSV.matrix.down[i,c(1,3)] <- -1
		RSV.matrix.down[i,2] <- 0
	}
	if (RSV.ids.down[i,2]==8){
		RSV.matrix.down[i,2:3] <- -1
		RSV.matrix.down[i,1] <- 0
	}
	if (RSV.ids.down[i,2]==9){
		RSV.matrix.down[i,1:3] <- -1
	}
}

RSV.matrix.updown <- rbind(RSV.matrix.up, RSV.matrix.down)


### Venn Diagrams ###

par(mfrow=c(1,2))
vennDiagram(RSV.matrix, names=c("RSVvLC", "RSVvZC", "ZCvLC"), circle.col=c("darkgreen", "orange", "blue"), main=paste(TISSUE, "Overlap in All DE Genes", sep=" "))
vennDiagram(RSV.matrix.updown, include=c("up", "down"), names=c("RSVvLC", "RSVvZC", "ZCvLC"), circle.col=c("darkgreen", "orange", "blue"), main=paste(TISSUE, "Overlap in DE Genes", sep=" "), sub="test")
quartz.save(paste("~/Desktop/", paste(TISSUE, "Venn RSV.pdf", sep=" "), sep=""), type="pdf", width=13, height=7)
