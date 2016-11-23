#############################################
## Rafi Pelossof, MSK 2014
#
# sexen - SEquential eXclusive ENrichment test
#

su.hypergeom <- function (binsize, ngenes, median.rank) {


	k <- median.rank
	nk <- ngenes - k
	p <- phyper(binsize/2, binsize, nk, k, lower.tail=FALSE)
	p[is.nan(p)] <- 1
	#if (k - nk > binsize) print('su.hypergeom ill conditioned')

	return (p)
}

##su.pvals <- function (ranks) {
##	binsize <- length(ranks[[1]])
##	ntotal <- length(unlist(ranks))
##	pvals <- lapply (ranks, function (rank) {su.hypergeom(binsize, ntotal, median(rank))})
##	return (unlist(pvals))
##}


su.stats <- function (x, y, binsize = 50) {
# vectors with:
# ordered genes lists by some metric (fdr), values(x) are lfc, names(x) are gene names

	genes <- intersect (rownames(x), rownames(y))
	x <- x[rownames(x) %in% genes, ]
	y <- y[rownames(y) %in% genes, ]


	y$rank <- as.numeric(1:nrow(y))
	print(head(x))
	print(head(y))

	ix <- seq(0, nrow(y), binsize)
	ranks <- list()
	pvals <- c()
	x.padj <- list()
	for (i in 1:(length(ix)-1)) {
		ranks[[i]] <- y[rownames(x)[(ix[i]+1):ix[i+1]], 'rank']
		x.padj[[i]] <- x[rownames(x)[(ix[i]+1):ix[i+1]], 'padj']
		pvals[i] <- su.hypergeom(binsize, nrow(x), median(ranks[[i]]))
	}
	return (list(ranks = ranks, x.padj = x.padj, pvals = pvals))
}


su.plot <- function (su.result, main = '') {
	ranks <- c(su.result$up$ranks, rev(su.result$down$ranks))
	ln.rank.up <- length(su.result$up$ranks)
	ln.rank.down <- length(su.result$down$ranks)
	pvals <- c(su.result$up$pvals, rev(su.result$down$pvals))
	ln.pval.up <- length(su.result$up$pvals)
	ln.pval.down <- length(su.result$down$pvals)
	x.padj <- c(sapply(su.result$up$x.padj, function(x) -log10(x), simplify = F), 
		rev(sapply(su.result$down$x.padj, function(x) -log10(x), simplify = F)))

	bin_size = length(su.result$up$ranks[[1]])

	par(mfrow = c(3, 1), mar = c(2, 2, 2, 2))
	boxplot (ranks, col = c(rep("#e41a1c", ln.rank.up), rep("#377eb8", ln.rank.down)), 
		pch = 16, font = 2, cex = 0.8, font.axis = 2, font.lab = 2, main = sprintf('%s, bin size %d', main, bin_size ))
	boxplot (x.padj, col = c(rep("#e41a1c", ln.rank.up), rep("#377eb8", ln.rank.down)), 
		pch = 16, font = 2, cex = 0.8, font.axis = 2, font.lab = 2,  main = sprintf('%s, padj', main ))
	abline(-log10(0.05), 0, "gray")
	#plot (-log10(pvals), col = c(rep("#e41a1c", ln.pval.up), rep("#377eb8", ln.pval.down)),
	#	pch = 16, font = 2, cex = 0.8, font.axis = 2, font.lab = 2)
	#abline(-log10(0.1), 0, col = 'gray')
	plot (pvals, col = c(rep("#e41a1c", ln.pval.up), rep("#377eb8", ln.pval.down)),
		pch = 16, font = 2, cex = 0.8, font.axis = 2, font.lab = 2,  main = sprintf('%s, pvals', main ))
	abline(0.05, 0, col = 'gray')
	#plot(x = rank.list.result$genes[ ,1], y = -log10(rank.list.result$genes[ ,2]), pch = 16,
	#	xlab = "log2FC, ylab = -log10(padj)", font.axis = 2, font.lab = 2)
	#abline(-log10(0.05), 0, col = "gray")
#	abline(length(genes)/2,0, lty=2, col = 'gray')
	

}


## main calling function
su.stats.updown <- function (x, y, binsize = 50) {
	# split list into two
	enrich.up <- su.stats (subset(x, lfc > 0) , subset(y, lfc > 0), binsize)
	enrich.down <- su.stats (subset(x, lfc < 0) , subset(y, lfc < 0), binsize)

	return (list (up = enrich.up, down = enrich.down))	
}


su.run <- function (x,y, padj.name = 'padj', lfc.name = 'log2FoldChange', binsize = 50 ) {
#expects two data frames. Rownames are genes, column names are set by default to padj and log2FoldChange
	#genes = intersect(rownames(x), rownames(y))

	x.su = x[ , c(padj.name, lfc.name)]
	y.su = y[ , c(padj.name, lfc.name)]
	colnames(x.su) = c('padj', 'lfc')
	colnames(y.su) = c('padj', 'lfc')
	x.su = x.su[order(x.su$padj), ]
	y.su = y.su[order(y.su$padj), ]

	stats = su.stats.updown(x.su, y.su, binsize = binsize)
}

su.limma <- function (tt.x, tt.y, binsize = 50) {
# test if tt.x is enriched in tt.y
	#clena tt.x
	dups = duplicated(tt.x$gene)
	tt.x <- tt.x[!dups,]
	dups = duplicated(tt.y$gene)
	tt.y <- tt.y[!dups,]
	x.ol = tt.x$logFC[order(tt.x$adj.P.Val)]
	names(x.ol) <- tt.x$gene[order(tt.x$adj.P.Val)]
	y.ol = tt.y$logFC[order(tt.y$adj.P.Val)]
	names(y.ol) <- tt.y$gene[order(tt.y$adj.P.Val)]
	ranks = su.rank.list(x.ol, y.ol, binsize = binsize)
	pvals = list(up = su.pvals(ranks[['up']]), down = su.pvals(ranks[['down']]))
	return (list(ranks = ranks, pvals = pvals))
}



su.example <- function () {
	#generate datasets by simulation
	n = 100;
	binsize = 20;
	x = data.frame (lfc = rnorm(n), fdr = runif(n), row.names = sprintf('g%d', 1:n))
	y = data.frame (lfc = x$lfc + rnorm(n)/5, fdr = x$fdr + runif(n)/2, row.names = rownames(x))


	y.control = rnorm(n)
	names(y.control) <- names(x)

	enrich.test = su.rank.list(x, y.test)
	enrich.control = su.rank.list(x, y.control)



	#load some existing dataset

}




