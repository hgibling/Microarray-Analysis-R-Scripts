##### Determine Overlapped Genes for Rosiglitazone #####

library(plyr)

ROSvLC.id <- cbind(ROSvLC.genes, 1)
ROSvZC.id <- cbind(ROSvZC.genes, 3)
ZCvLC.id <- cbind(ZCvLC.genes, 5)

ROS.genes <- as.data.frame(rbind(ROSvLC.id, ROSvZC.id, ZCvLC.id))

colnames(ROS.genes) <- c("GeneID", "Pair")

ROS.ids <- ddply(ROS.genes, "GeneID", summarize, Pair=sum(Pair))

ROS.LC <- which(ROS.ids[,2]==1)
ROS.ZC <- which(ROS.ids[,2]==3)
ZC.LC <- which(ROS.ids[,2]==5)
ROS.LC.ZC <- which(ROS.ids[,2]==4)
ROS.LC.ZCLC <- which(ROS.ids[,2]==6)
ROS.ZC.ZCLC <- which(ROS.ids[,2]==8)
ROS.all <- which(ROS.ids[,2]==9)


##### Venn Diagrams #####

library(VennDiagram)

draw.triple.venn(area1=length(ROSvLC.genes), area2=length(ROSvZC.genes), area3=length(ZCvLC.genes), n12=length(ROS.LC.ZC), n23=length(ROS.ZC.ZCLC), n13=length(ROS.LC.ZCLC), n123=length(ROS.all), category=c("ROSvLC", "RSVvLC", "ZCvLC"), scaled=FALSE)
#doesn't work with large overlap

library("venneuler")

venn.ROS <- venneuler(c(ROSvLC=length(ROS.LC), ROSvZC=length(ROS.ZC), ZCvLC=length(ZC.LC), "ROSvLC&ROSvZC"=length(ROS.LC.ZC), "ROSvLC&ZCvLC"=length(ROS.LC.ZCLC), "ROSvZC&ZCvLC"=length(ROS.ZC.ZCLC), "ROSvLC&ROSvZC&ZCvLC"=length(LC.all)))
plot(venn.ROS)
quartz.save("~/Desktop/Blood Venn ROS Proportional.pdf", type="pdf")
#works with large overlap, but can't unscale

library(limma)

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

ROS.counts <- vennCounts(ROS.matrix)

vennDiagram(ROS.counts, names=c("ROSvLC", "ROSvZC", "ZCvLC"), circle.col=c("darkgreen", "orange", "blue"), main="Overlap in DE Genes")
quartz.save("~/Desktop/Blood Venn ROS.pdf", type="pdf")

m <- matrix(nrow=5, ncol=3, c(0,1,0,-1,1,0,0,-1,-1,0,1,1,-1,0,1))

mcount<- vennCounts(m, include="both")
vennDiagram(mcount, include="both")


library("Vennerable")



##### Determine Overlapped Genes for Resveratrol #####

library(plyr)

RSVvLC.id <- cbind(RSVvLC.genes, 1)
RSVvZC.id <- cbind(RSVvZC.genes, 3)
ZCvLC.id <- cbind(ZCvLC.genes, 5)

RSV.genes <- as.data.frame(rbind(RSVvLC.id, RSVvZC.id, ZCvLC.id))

colnames(RSV.genes) <- c("GeneID", "Pair")

RSV.ids <- ddply(RSV.genes, "GeneID", summarize, Pair=sum(Pair))

RSV.LC <- which(RSV.ids[,2]==1)
RSV.ZC <- which(RSV.ids[,2]==3)
ZC.LC <- which(RSV.ids[,2]==5)
RSV.LC.ZC <- which(RSV.ids[,2]==4)
RSV.LC.ZCLC <- which(RSV.ids[,2]==6)
RSV.ZC.ZCLC <- which(RSV.ids[,2]==8)
RSV.all <- which(RSV.ids[,2]==9)


##### Venn Diagrams #####

library(VennDiagram)

draw.triple.venn(area1=length(RSVvLC.genes), area2=length(RSVvZC.genes), area3=length(ZCvLC.genes), n12=length(RSV.LC.ZC), n23=length(RSV.ZC.ZCLC), n13=length(RSV.LC.ZCLC), n123=length(RSV.all), category=c("ROSvLC", "RSVvLC", "ZCvLC"), scaled=FALSE)
#doesn't work with large overlap

library("venneuler")

venn.RSV <- venneuler(c(RSVvLC=length(RSV.LC), RSVvZC=length(RSV.ZC), ZCvLC=length(ZC.LC), "RSVvLC&RSVvZC"=length(RSV.LC.ZC), "RSVvLC&ZCvLC"=length(RSV.LC.ZCLC), "RSVvZC&ZCvLC"=length(RSV.ZC.ZCLC), "RSVvLC&RSVvZC&ZCvLC"=length(RSV.all)))
plot(venn.RSV)
quartz.save("~/Desktop/Blood Venn RSV Proportional.pdf", type="pdf")
#works with large overlap, but can't unscale

library(limma)

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

RSV.counts <- vennCounts(RSV.matrix)

vennDiagram(RSV.counts, names=c("RSVvLC", "RSVvZC", "ZCvLC"), circle.col=c("red", "purple", "blue"), main="Overlap in DE Genes")
quartz.save("~/Desktop/Blood Venn RSV.pdf", type="pdf")

m <- matrix(nrow=5, ncol=3, c(0,1,0,-1,1,0,0,-1,-1,0,1,1,-1,0,1))

mcount<- vennCounts(m, include="both")
vennDiagram(mcount, include="both")


library("Vennerable")

