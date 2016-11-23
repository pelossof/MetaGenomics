
#Exen - exclusive enrichment of one list in another

exen.hypergeom <- function (binsize, ngenes, median.y) {

	k <- median.y
	n <- ngenes - k
	m <- binsize
	x <- binsize/2
	p <- phyper(x, m, n, k, lower.tail=FALSE)

	return (p)
}

exen.rank.en <- function (x, y, binsize = 100) {
# vectors with:
# ordered genes lists by some metric (fdr), values(x) are lfc, names(x) are gene names
	genes <- intersect (names(x), names(y))

	x <- x[genes]
	y <- y[genes]

	y[,'rank'] <- 1:length(y)

	ix <- seq(1, length(genes), binsize)

	b = list()
	for (i in 1:(length(ix)-1)) {					## all the subsequent boxes
       b[[i]]$ranks = y[rownames(group2)[ix[i]:(ix[i+1]-1)],'rank']
       b[[i]]$pval = 
	}
	return (b)
}

	

}

exen.rank.list <- function (x, y) {
	enrich.up <- exen.rank.en (x[x > 0] , y[y > 0])
	enrich.down <- exen.rank.en (x[x < 0] , y[y < 0])

}

exen.split <- function (x, y) {

}


exen.deseq <- function (x, x.fc, y, y.fc) {
# takes two results of DESeq and checks until where is x enriched in y


}




plotBoxEnrich <- function(group1, group1name="Group1", group2, group2name="Group2", binsize = 100, date = "test"){
## e.g. group1 - TCGA (rank column based on padj), group2 - Timing
#plotBoxEnrich(res_tc19, "TCGA", timing$mir$de, "Timing", 20, "150525")

binsize = 100;
	group1 <- removeNA(group1)  
	group2 <- removeNA(group2)

	genes <- intersect(rownames(group2), rownames(group1))
	group1 <- group1[genes, ]
	group2 <- group2[genes, ]
	group1 <- group1[order(group1$padj), ]
	group2 <- group2[order(group2$padj), ]

	group1$rank = 1:nrow(group1)

	selectlist <- subset(group2, group2$padj < 0.2) # & abs(group2$log2FoldChange) > 1)

	ix = seq(1, length(genes), binsize)

	b = list()
	b[[1]] <- group1[rownames(selectlist), 'rank']  	## the selected feature list
	
	for (i in 1:(length(ix)-1)) {					## all the subsequent boxes
       b[[i+1]] = group1[rownames(group2)[ix[i]:(ix[i+1]-1)],'rank']
	}

	colors <- c("red", rep("white", length(b)-1))
	labels <- c("Signature", paste("Bin ", (2:length(b))-1, sep=""))

	pdf(paste("results/", date, "_boxEnrich_", group2name, "_", group1name, "_binsize", binsize, ".pdf", sep=""))
		#par(mar=c(7.1,4.1,3.1,3.1))  # (bottom, left, top, right)

		boxplot(b, col=colors, names = labels, las=2, outline=FALSE, ylab="TCGA rank", 
			main=paste("Enrichment analysis of ", length(genes), 
				" shared features\nselected from TCGA and Timing. Binsize: ", binsize, ".", sep=""))

		abline(length(genes)/2,0, lty=2, col = 'gray')
	dev.off()
}

