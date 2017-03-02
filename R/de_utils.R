
library(DESeq2)
library (GSA)



pval2sym = function (pvals) {
	syms = rep(NA, length(pvals))
	syms[pvals > 0.1] = "ns"
	syms[pvals <= 0.1] = "-"
	syms[pvals <= 0.05] = "*"
	syms[pvals <= 0.01] = "**"
	syms[pvals <= 0.001] = "***"
	syms[pvals <= 0.0001] = "****"
	names(syms) = names(pvals)
	return (syms)
}


pca.plot = function(dat, pc.ix1 = 1, pc.ix2  = 2, ann = colnames(dat), ann.cex = 0.5, ...) {
	#columns are samples, rows are dimensions
	pca = prcomp (t(dat))
	pred = predict(pca, t(dat))
	plot(pred[, pc.ix1], pred[,pc.ix2], ...)
	text(pred[, pc.ix1], pred[,pc.ix2], ann, cex = ann.cex, pos = 3)
}


########################
##  GSEA code         ##
########################
write.gsea = function(genes, filename, name = "") {
# for gene lists use extension .gmx
	write.table(c(name, unname(unlist(genes))), file = filename, row.names = F, col.names = F, quote = F)
}

write.gsea.de = function (de, filename, lfc = 'log2FoldChange', padj = 'padj') {
# for ranked lists use extension .rnk
# de is the resulting table of DESeq analysis: $res	
	de = de[!is.na(de[,padj]), ]
	de$score = sign(de[,lfc]) * -log10(de[,padj])
	de = subset(de, !is.na(score))
	#de$score[is.na(de$score)] = 0
	score = sprintf('%.8f', de$score)
	names(score) = rownames(de)
	write.table(score, file = filename, row.names = T, col.names = F, quote = F, sep = '\t')
}

run.gsea = function(rnk.fname, sig.fname, out.path, analysis.name, set.max = 500, set.min = 15, nperm = 1000, collapse = 'false', mem = '2g', gsea.jar = '/Users/pelossof/tools/gsea2-2.2.1.jar', wait = F) {
	cmd = list()
	cmd$cmd = sprintf('java -Xmx%s -cp %s xtools.gsea.GseaPreranked', mem, gsea.jar)
	cmd$rpt = sprintf('-rpt_label %s', analysis.name)
	cmd$gmx = sprintf('-gmx %s', sig.fname) 
	#cmd$chip = sprintf('-chip %s', paths.data('gsea/chip_files/GENE_SYMBOL.chip'))
	cmd$out = sprintf('-out %s', out.path)
	cmd$collaps = sprintf('-collapse %s', tolower(as.character(collapse))) 
	cmd$mode = '-mode Max_probe'
	cmd$norm = '-norm meandiv'
	cmd$nperm = sprintf('-nperm %d', nperm)
	cmd$rnk = sprintf('-rnk %s', rnk.fname)
	cmd$scoring = '-scoring_scheme weighted' 
	#cmd$help = '-help false' 
	cmd$sym = '-include_only_symbols true' 
	cmd$sets = '-make_sets true' 
	cmd$plot = '-plot_top_x 20' 
	cmd$seed = '-rnd_seed timestamp' 
	cmd$max = sprintf('-set_max %d', set.max) 
	cmd$min = sprintf('-set_min %d', set.min) 
	cmd$zip = '-zip_report false' 
	cmd$gui = '-gui false'

	cmd.run = paste(cmd, sep = '', collapse = ' ')
	print(sprintf('Running GSEA:> %s', cmd.run))
	system(cmd.run, wait)
}



read.gmt.all = function() {
	gmts = list.files(paths.data('gsea/msigdb_v5.0_GMTs'), pattern = '.*all.*symbol.*', full.names = T)
	gsea.lists = lapply (gmts, read.gmt1)
	names(gsea.lists) = basename(gmts)
	gsea.all = unlist(gsea.lists, recursive = F)
	return(gsea.all)
}

run.gsea.all = function (rnk.fname, gsea.gmt.path = paths.data('gsea/msigdb_v5.0_GMTs'), pattern = '.*all.*symbol.*', mc.cores = 8) {
	gmts = list.files(gsea.gmt.path, pattern = pattern, full.names = T)
	mclapply (gmts, function (gmt) {run.gsea(rnk.fname, sig.fname = gmt, out.path = gsub('.rnk','', rnk.fname, fixed = T), analysis.name = basename(gmt))}, mc.cores = mc.cores)
}


read.gmt1 = function (filename) {
	x = GSA.read.gmt(filename)
	names(x$genesets) = x$geneset.names
	return(x$genesets)
}

read.gmt = function (filename) {
	x = read.delim(filename, header = F, stringsAsFactors = F, row.names = 1)
	y = lapply(1:nrow(x), function (r) as.character(x[r, 2:ncol(x)]))
	names(y) = rownames(x)
	return (y)
}

#### code for parsing GSEA output files - need to standardize
html2tab = function (url) { 
	doc = htmlParse(url)
	tableNodes = getNodeSet(doc, "//table")
	tb = readHTMLTable(tableNodes[[1]],header=F, as.data.frame = F)
	tab = as.data.frame(Reduce(cbind, tb), stringsAsFactors = F)
	colnames(tab) = tb.colnames
	tab.colnames = c('', 'sig', 'details', 'size','es','nes', 'pval', 'fdr','fwer', 'rank_at_max', 'leading_edge')
	colnames(tab) = tab.colnames
	tab = tab[,-1]

	tab = transform(tab,  sig = as.character(sig), 
						details = as.character(details), 
						size = as.numeric(size),
						es = as.numeric(es), 
						nes = as.numeric(nes),
						pval = as.numeric(pval), 
						fdr = as.numeric(fdr), 
						fwer = as.numeric(fwer), 
						rank_at_max = as.numeric(rank_at_max), 
						leading_edge = as.character(leading_edge))
	return (tab)
}

loadGSEAtable = function(folder) {
print(folder)
	neg.file = list.files(pattern = 'gsea_report.*neg.*.html', path = folder, full.names = T)
	pos.file = list.files(pattern = 'gsea_report.*pos.*.html', path = folder, full.names = T)
	
	return(list(pos = html2tab(pos.file), neg = html2tab(neg.file)))
}
#folders = list.files(pattern = '*symbols.gmt.GseaPreranked*')
#results = lapply(folders, loadGSEAtable)

#### END: code for parsing GSEA files - need to standardize



plotSignature = function (expr, ds, signatures.list, x, ann) {
	sig.names = names(signatures.list)
	for(sig.name in sig.names) {
		print(sprintf('plotting %s', sig.name))
		sig.genes = intersect(rownames(expr), signatures.list[[sig.name]])
		print(sig.genes)
		if(length(sig.genes) < 2 )
			next
		expr.g = expr[sig.genes, ]
		o = order(colMeans(expr.g, na.rm = T))
		expr.o = expr.g[, o]
		x.o = x[o]
		ann.o = ann[[1]][o]
		m = colMeans(expr.g, na.rm = T)
		m.o = m[o]
		plot(x, m, main = sig.name)
		f  = lm(m ~ x)
		abline(f)
		sf = summary(f)
		legend('topleft', legend = strsplit(sprintf('pvalue=%.5f,adj.r2=%.5f,coeff=%.5f',  sf$coefficients[2, 4], sf$adj.r.squared, sf$coefficients[2, 1]), ',')[[1]] )
		boxplot(m.o ~ ann.o, main = sig.name)
		for (l in 1:length(levels(ann.o)))
			points(jitter(rep(l, sum(ann.o == levels(ann.o)[l])), factor = 1, amount = 0.5/length(levels(ann.o))), m.o[ann.o == levels(ann.o)[l]], pch=19, col = rgb(0,0,1,0.5))

print('running aheatmap')
expr.o[is.na(expr.o)] = 0
		aheatmap(expr.o, annCol = list(var1 = ann.o), scale = 'row', Colv=NA, main = sig.name)

print('running volcano')
#		volcano.plot(ds$res, sig.genes, lfc.name = 'log2FoldChange', padj.name = 'padj', sig.fdr = 0.2, xlim = c(-2,2),ylim = c(0,10), main=sprintf('yng vs old\n%s', sig.name), arrows=T, significant.only=F)
		volcano.plot(ds$res, sig.genes, lfc.name = 'log2FoldChange', padj.name = 'padj', sig.fdr = 0.2, xlim = c(-2,2),ylim = c(0,10), main=sprintf('%s', sig.name), arrows=T, significant.only=F)
print('done volcano')

	}
}



oncoprinter.pats = function (pats, genes, cna, maf, filename, maf.pat.col = "Tumor_Sample_Barcode", maf.gene.col= "Hugo_Symbol") { 

	cn.g = rownames(cna)
	cn.p = colnames(cna)
	maf.p = maf[,maf.pat.col]

	pats = Reduce(intersect, list(pats, cn.p, maf.p))  ## patients that have all the data

	out = c()
	for (p in pats) {
		print(sprintf('processing %s', p ))
		mafrows.pat = maf[maf.p %in% p, ]
		cnas.pat = cna[, p]
		for (g in genes) {
		#	print(sprintf('%s/%s',p,g))
			if(g %in% cn.g) {
				cnas = cnas.pat[cn.g == g]
				if (cnas == 1) {
					out = rbind (out , c(p, g, 'AMP', 'CNA'))
				} else if (cnas == 2) {
					out = rbind (out , c(p, g, 'GAIN', 'CNA'))
				} else if (cnas == -1) {
					out = rbind (out, c(p, g, 'HETLOSS','CNA'))
				} else if (cnas == -2) {
					out = rbind (out, c(p, g, 'HOMDEL','CNA' ))
				}
			}

			nmut = sum(mafrows.pat[, maf.gene.col] == g)
			if (nmut == 1) {
				out = rbind (out, c(p, g, 'A1A', 'MISSENSE' ))
			} else if (nmut == 2) {
				out = rbind (out, c(p, g, 'A1A', 'TRUNC' ))
			}
		}
	}
	print(out)
	colnames(out) = c('Sample','Gene','Alteration','Type')
	write.table(out, filename, quote = F, row.names = F, col.names = T, sep='\t')
}




##
## create a plot for OncoPrinter
##
oncoprinter = function(dat, tags = NA, filename = 'oncoprinter.txt') {
	#dat: patients x genes
	#	matrix entries:
	#	  0: skip
	#	  1: 1 - block red (amp)
	#	  2: 2 - block blue (homdel)
	#	  3: 
	if (is.na(tags[1]))
		tags = rep('AMP', ncol(dat))

	names(tags) = rownames(dat)

	out = c()
	for (i in 1:ncol(dat)) {
		cn = colnames(dat)[i]
		nm = rownames(dat)[dat[,i]==1]
		nm = nm[!is.na(nm)]
		out = rbind(out, cbind(nm, rep(cn, length(nm)), tags[nm]))
	}
	## pCR will be in RED
	#out[out[,2]=='pCR', 3] = 'AMP'
	colnames(out) = c('Sample','Gene','Alteration')
	write.table(out, filename, quote = F, row.names = F, col.names = T, sep='\t')
}

##
##  DEseq2 code
##
deseq2 <- function (ct, factors,...) {
	colData <- data.frame(condition=factors, row.names=colnames(ct))
	dds <- DESeqDataSetFromMatrix (countData = ct, colData = colData, design = ~ condition)
	colData(dds)$condition <- factor(colData(dds)$condition, levels = levels(factors))
	dds <- DESeq(dds,...)
	res <- results(dds)
	o <- order(res$padj)
	res <- res[o,]

	return (list(res = res, dds = dds))
}

deseq2plots2 = function (ct, is.test, study, col.null, col.test, name.null, name.test, volcano.signatures = list(), v , run_gsea = F, basemean.min = 50, ...) {
	return (deseq2plots(ct, !is.test, study, col.null, col.test, name.null, name.test, volcano.signatures, v , run_gsea, basemean.min, ...))
}

deseq2plots = function (ct, is.null, study, col.null, col.test, name.null, name.test, volcano.signatures = list(), v , run_gsea = F, basemean.min = 50, ...) {
	factors.all <- factor( is.null, levels = c(TRUE, FALSE), labels = c('null', 'test') )
	ds.test <- deseq2(round(ct), factors.all)
	ds.test$res = res2gene(ds.test$res)
	plotDEseq2Results(ds.test$dds, ds.test$res, ver.path(sprintf('%s_rnaseq_dispersion.pdf', study),v))
	write.csv(ds.test$res, ver.path(sprintf('%s_de.csv', study), v))
	pdf(ver.path(sprintf('%s_pca_volcano.pdf',study), v))
	pca.plot(log2(1+ct), col = ifelse(is.null, col.null, col.test), pch = 19, ann.cex = 0.5)
	legend('topleft', c(name.test, name.null), col = c(col.test, col.null), pch = 19)
	for(vol.sig in names(volcano.signatures)) {
		volcano.plot(ds.test$res, volcano.signatures[[vol.sig]], main = sprintf('%s, %s', study, vol.sig), ylim = c(0, max(-log10(ds.test$res$padj + 1e-100), na.rm = T)), ...)
		legend('top', sprintf('<-- %s high | %s high -->', name.null, name.test))
	}
	dev.off()
	if (run_gsea) {		
		filename.rnk = ver.path(sprintf('%s.rnk', study), v)
		write.gsea.de(subset(ds.test$res, baseMean > basemean.min), filename.rnk, lfc = 'log2FoldChange', padj = 'padj')
		run.gsea.all(filename.rnk)
		run.gsea(filename.rnk, sig.fname = paths.data('gsea/signatures/sig_stroma_20160914.gmt'), out.path = gsub('.rnk','', filename.rnk, fixed = T), analysis.name = 'sig_stroma')
	}
	return (ds.test)
}


diffExpr <- function (arr, lab, arr.type, intercept = TRUE, p.value = 0.1, number = 1000) {
#arr is a probe intensity matrix where each row is a probe, and each column is a sample
#lab is the label to check diff expression for

	if (intercept) {
		#print("fitting a model with an intercept")
		design <- cbind(intercept = rep(1, length(lab)), lab = lab)
		fit <- lmFit(arr, design)
		efit <- eBayes(fit)
		tt <- topTable(efit, p.value = p.value, number = number, coef="lab", adjust="fdr")
	} else {
		#print("fitting a model without an intercept")
		fit <- lmFit(arr, lab)
		efit <- eBayes(fit)
		tt <- topTable(efit, p.value = p.value, number = number, coef="lab", adjust="fdr")
	}

#	if (nrow(tt)>0)
#		tt$gene <- probes.to.symbols(tt$ID, arr.type)

	return (tt)
}


plotDEseq2Results <- function (dds, res, basename) {
	#pdf(sprintf('results/deseq2_%s.pdf', basename)); 
	pdf(basename)
	plotDispEsts(dds);
	DESeq2::plotMA(dds);
	hist(res$pvalue, breaks=100, col="skyblue", border="slateblue", main=""); 
	dev.off();
}

plotReplicatesLFC <- function (res1, res2, padj_thresh, basename) {
	pdf(basename)
	genes1 <- subset(res1, padj < padj_thresh)
	genes2 <- subset(res2, padj < padj_thresh)

	genes <- intersect(rownames(genes1), rownames(genes2))

	print(sprintf('found %d genes in the intersection', length(genes)))

	plot(res1[genes, 'log2FoldChange'], res2[genes, 'log2FoldChange'], xlim=c(-3,3), ylim=c(-3,3))
	abline(a=0, b=1);
	dev.off()
}

listEnrichment <- function (list1, list2, bin.size=200) {

}



limma2deseq <- function (res) {
	n <- colnames(res)
	n <- gsub('logFC', 'log2FoldChange', n)
	n <- gsub('adj.P.Val', 'padj', n)
	n <- gsub('AveExpr', 'baseMean', n)
	colnames(res) <- n
	return (res)
}

removeNA <- function (res) {
	return (res[!is.na(res$padj), ])
}

gene.rownames <- function (de.list, gene.names) {
	de.list = de.list[!duplicated(toupper(de.list[,gene.names])) &! is.na(de.list[,gene.names]),]
	rownames(de.list) = toupper(de.list[,gene.names])
	return (de.list)
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



volcano.plot.old <- function (x, genes, main = '', lfc.name = 'logFC', padj.name = 'adj.P.Val', 
		sig.fdr = 0.05, x.offset = 0.01*(max(xlim)-min(xlim)), y.offset = 0.025*(max(ylim)-min(ylim)), xlim = c(-3,3), ylim = c(0,1.8), 
		legend.text = c("Significant", "Not Significant", "All genes")) {
	padj = padj.name
	lfc = lfc.name

	bold = 2

	x <- x[!is.na(x[,lfc]) & !is.na(x[,padj]), ]
	ix = rownames(x) %in% genes
	sig.ix = ix & (x[,padj] < sig.fdr)
	plot(x[, lfc], -log10(x[, padj]), col = "gray", main = main, cex = 0.5, pch = 20, xlim = xlim, ylim = ylim)
	if (sum(ix) > 0) {
		points(x[ix, lfc], -log10(x[ix, padj]), col = "dodgerblue", pch = 20)
		text(x[ix, lfc]+x.offset, -log10(x[ix, padj])+y.offset, rownames(x[ix, ]), cex = 0.8, font = bold )
	}
	if (sum(sig.ix) > 0) {
		points(x[sig.ix, lfc], -log10(x[sig.ix, padj]), col = 'red', pch = 20)
#		text(x[sig.ix, lfc]+x.offset, -log10(x[sig.ix, padj])+y.offset, rownames(x[sig.ix, ]), cex = 0.5, font = bold )
	}
	abline(-log10(sig.fdr), 0)
#	text(min(xlim) + 0.01*(max(xlim)-min(xlim)), -log10(sig.fdr)+y.offset,sprintf('fdr < %.3f', sig.fdr), col='red', cex = 0.5, font = bold)
	text(min(xlim) + 0.01*(max(xlim)-min(xlim)), -log10(sig.fdr)+y.offset,sprintf('fdr < %.3f', sig.fdr), col='red', cex = 0.5, font = bold)

	legend("bottomleft", legend.text, col = c("red", "dodgerblue", "gray"), pch = c(20,20,20), cex = 0.5)
}

volcano.plot <- function (x, genes, main = '', lfc.name = 'log2FoldChange', padj.name = 'padj', 
		sig.fdr = 0.05, x.offset = 0.01*(max(xlim)-min(xlim)), y.offset = 0.025*(max(ylim)-min(ylim)), xlim = c(-3,3), ylim = c(0,1.8), 
		fx = function(x){x}, fy=function(y){-log10(y)},
		xlab = "log2FoldChange", ylab = "-log10(padj)",
		legend.text = c(paste("Significant fdr < ", sig.fdr), "Not Significant", "All genes"), arrows = FALSE, max.genes = 50, 
		significant.only = T,
		point.colors = c(), 
		samp.percent = 1, ...) {
	padj = padj.name
	lfc = lfc.name

	bold = 2

	if(length(point.colors) > 0)
		names(point.colors) = genes
	x <- x[!is.na(x[,lfc]) & !is.na(x[,padj]), ]
	x[,lfc] <- fx(x[,lfc])
	x[,padj] <- fy(x[,padj])
	x <- x[order(x[,padj], decreasing=T),]
	ix = rownames(x) %in% genes
	x.g = x[ix, ]
	sig.ix = x.g[,padj] > fy(sig.fdr)
	ix.rnd = sample(nrow(x), round(nrow(x)*samp.percent), replace=F)
	plot(x[ix.rnd, lfc], x[ix.rnd, padj], xlab=xlab, ylab = ylab, col = "gray", main = main, pch = 20, xlim = xlim, ylim = ylim, cex = 1)
	if (sum(ix) > 0) {
		if(length(point.colors) > 0)
			points(x.g[, lfc], x.g[, padj], col = point.colors[rownames(x.g)], ...)
		else
			points(x.g[, lfc], x.g[, padj], col = "dodgerblue", ...)


		if (arrows==TRUE) {
			#split to up regulated and down regulated
			if (significant.only) {
				up = x.g[(x.g[, lfc]>0) & sig.ix, ]
				down = x.g[(x.g[, lfc]<0) & sig.ix, ]
			} else {
				up = x.g[x.g[, lfc]>0, ]
				down = x.g[x.g[, lfc]<0, ]
			}
			dy=(max(ylim)-min(ylim))/65
			if ( nrow(up) > 0) {
				for (i in 1:min(max.genes, nrow(up))) {

					text(max(xlim)-0.05*(max(xlim)-min(xlim)), max(ylim)-dy*(i-1), rownames(up)[i], cex=0.5, font=1, adj=c(0,0.5), col=ifelse(up[i, padj] > fy(sig.fdr), 'red','dodgerblue'))
					segments(max(xlim)-0.07*(max(xlim)-min(xlim)), max(ylim)-dy*(i-1), up[i,lfc], up[i, padj], lwd=0.5, col=rgb(0,0,0,0.2))
					segments(max(xlim)-0.06*(max(xlim)-min(xlim)), max(ylim)-dy*(i-1), 
						max(xlim)-0.07*(max(xlim)-min(xlim)), max(ylim)-dy*(i-1), lwd=0.5, col=rgb(0,0,0,0.2))
				}
			}
			if ( nrow(down) >0 ) {
				for (i in 1:min(max.genes, nrow(down))) {
					text    (min(xlim)+0.05*(max(xlim)-min(xlim)), max(ylim)-dy*(i-1), rownames(down)[i], cex=0.5, font=1, adj=c(1,0.5), col=ifelse(down[i, padj] > fy(sig.fdr),  'red', 'dodgerblue'))
					segments(min(xlim)+0.07*(max(xlim)-min(xlim)), max(ylim)-dy*(i-1), down[i, lfc], down[i, padj], lwd=0.5, col=rgb(0,0,0,0.2))
					segments(min(xlim)+0.06*(max(xlim)-min(xlim)), max(ylim)-dy*(i-1), min(xlim)+0.07*(max(xlim)-min(xlim)), max(ylim)-dy*(i-1), lwd=0.5, col=rgb(0,0,0,0.2))
				}			
			}
		} #else {
		#	text(x[ix, lfc]+x.offset, x[ix, padj]+y.offset, rownames(x[ix, ]), cex = 0.8, font = bold )
		#}
	}
	if (sum(sig.ix) > 0 & length(point.colors) == 0) {
		points(x.g[sig.ix, lfc], x.g[sig.ix, padj], col = 'red', pch = 20)
	#	text(x[sig.ix, lfc]+x.offset, -log10(x[sig.ix, padj])+y.offset, rownames(x[sig.ix, ]), cex = 0.5, font = bold )
	}
	#			
#	abline(-log10(sig.fdr), 0)
#	text(min(xlim) + 0.01*(max(xlim)-min(xlim)), -log10(sig.fdr)+y.offset,sprintf('fdr < %.3f', sig.fdr), col='red', cex = 0.5, font = bold)
#	text(min(xlim) + 0.01*(max(xlim)-min(xlim)), -log10(sig.fdr)+y.offset,sprintf('fdr < %.3f', sig.fdr), col='red', cex = 0.5, font = bold)


#	legend("bottomleft", legend.text, col = c("red", "dodgerblue", "gray"), pch = c(20,20,20), cex = 0.7)
}



ma.plot <- function (x, genes, main = '', lfc.name = 'logFC', padj.name = 'adj.P.Val', 
		sig.fdr = 0.05, x.offset = 0.01*(max(xlim)-min(xlim)), y.offset = 0.025*(max(ylim)-min(ylim)), xlim = c(-3,3), ylim = c(0,1.8), 
		fx = function(x){x}, fy=function(y){-log10(y)},
		xlab = "log2FoldChange", ylab = "-log10(padj)",
		legend.text = c(paste("Significant fdr < ", sig.fdr), "Not Significant", "All genes"), arrows = FALSE, max.genes = 50, 
		significant.only = T, ...) {
	padj = padj.name
	lfc = lfc.name

	bold = 2

	x <- x[!is.na(x[,lfc]) & !is.na(x[,padj]), ]
	x[,lfc] <- fx(x[,lfc])
	x[,padj] <- fy(x[,padj])
	x <- x[order(x[,padj], decreasing=T),]
	ix =  rownames(x) %in% genes
	x.g = x[ix, ]
	sig.ix = x.g[, 'padj'] < sig.fdr
	plot(x[, lfc], x[, padj], xlab=xlab, ylab = ylab, col = "gray", main = main, cex = 0.5, pch = 20, xlim = xlim, ylim = ylim,...)
	if (sum(ix) > 0) {
		points(x.g[, lfc], x.g[, padj], col = "dodgerblue", pch = 20)

		if (arrows==TRUE) {
			#split to up regulated and down regulated
			if (significant.only) {
				up = x.g[(x.g[, lfc]>0) & sig.ix, ]
				down = x.g[(x.g[, lfc]<0) & sig.ix, ]
			} else {
				up = x.g[x.g[, lfc]>0, ]
				down = x.g[x.g[, lfc]<0, ]
			}
			dy=(max(ylim)-min(ylim))/65
			if ( nrow(up) ) {
				#print(up)
				for (i in 1:min(max.genes, nrow(up))) {
					text(max(xlim)-0.05*(max(xlim)-min(xlim)), max(ylim)-dy*(i-1), rownames(up)[i], cex=0.5, font=1, adj=c(0,0.5), col=ifelse(up[i, 'padj'] < sig.fdr, 'red','dodgerblue'))
					segments(max(xlim)-0.07*(max(xlim)-min(xlim)), max(ylim)-dy*(i-1), up[i,lfc], up[i, padj], lwd=0.5, col=rgb(0,0,0,0.2))
					segments(max(xlim)-0.06*(max(xlim)-min(xlim)), max(ylim)-dy*(i-1), 
						max(xlim)-0.07*(max(xlim)-min(xlim)), max(ylim)-dy*(i-1), lwd=0.5, col=rgb(0,0,0,0.2))
				}
			}
			if ( nrow(down) ) {
				#print(down)
				for (i in 1:min(max.genes, nrow(down))) {
					text    (min(xlim)+0.05*(max(xlim)-min(xlim)), max(ylim)-dy*(i-1), rownames(down)[i], cex=0.5, font=1, adj=c(1,0.5), col=ifelse(down[i, 'padj'] < sig.fdr,  'red', 'dodgerblue'))
					segments(min(xlim)+0.07*(max(xlim)-min(xlim)), max(ylim)-dy*(i-1), down[i, lfc], down[i, padj], lwd=0.5, col=rgb(0,0,0,0.2))
					segments(min(xlim)+0.06*(max(xlim)-min(xlim)), max(ylim)-dy*(i-1), min(xlim)+0.07*(max(xlim)-min(xlim)), max(ylim)-dy*(i-1), lwd=0.5, col=rgb(0,0,0,0.2))
				}			
			}
		} else {
			text(x[ix, lfc]+x.offset, x[ix, padj]+y.offset, rownames(x[ix, ]), cex = 0.8, font = bold )
		}
	}
	if (sum(sig.ix) > 0) {
		points(x.g[sig.ix, lfc], x.g[sig.ix, padj], col = 'red', pch = 20)
#		text(x[sig.ix, lfc]+x.offset, -log10(x[sig.ix, padj])+y.offset, rownames(x[sig.ix, ]), cex = 0.5, font = bold )
	}
#	abline(-log10(sig.fdr), 0)
#	text(min(xlim) + 0.01*(max(xlim)-min(xlim)), -log10(sig.fdr)+y.offset,sprintf('fdr < %.3f', sig.fdr), col='red', cex = 0.5, font = bold)
#	text(min(xlim) + 0.01*(max(xlim)-min(xlim)), -log10(sig.fdr)+y.offset,sprintf('fdr < %.3f', sig.fdr), col='red', cex = 0.5, font = bold)


	legend("bottomleft", legend.text, col = c("red", "dodgerblue", "gray"), pch = c(20,20,20), cex = 0.7)
}




#' @export
plot.ecdf <- function (x ,y, xlim=c(min(x, na.rm=T), max(x, na.rm=T)), ylim = c(0,1), title='ECDF', xlab='Values', ylab='Fraction', legends=NA,...) {
 x = x[!is.na(x)]
 y = y[!is.na(y)]

 p.value.over <- ks.test (x, y, alternative="less")$p.value
 p.value.under <- ks.test (x, y, alternative="greater")$p.value

 ylim = c(0, 1)
 plot (0, 0, xlim=xlim, ylim = ylim,  type="n", xlab=xlab, ylab=ylab, main=title,...)
 plot (stepfun (sort (x), (0:length(x))/length(x)), col="red", add=TRUE, pch='')
 plot (stepfun (sort (y), (0:length(y))/length(y)), col="black", add=TRUE, pch='')

 if (hasArg(legends))
   legend ("bottomright", legend=legends, col=c("red", "black"), lty=1, bty="n")
 legend ("topleft", legend=sprintf ("KS test \nOver: %1.2e\nUnder: %1.2e", p.value.over, p.value.under), lty=0, bty="n")

}

#' @export
add.to.ecdf <- function (data, col="blue", ...) {
 plot (stepfun (sort (data), (0:length (data))/length (data)), col=col, add=TRUE, pch='', ...)
}

boxplot.jitter = function(x, y, jit = 0.4, pch = 20, col = rgb(0,0,1,.5), ...) {
	df = as.data.frame(cbind(x = x, y = y, jit = 1+unlist(jitter(as.numeric(x), jit) )))
	#plot.ecdf(unlist(subset(df, x == 1, y)), unlist(subset(df, x == 0, y)), ...)
	boxplot(y ~ x, data = df, outline = F,...)
	points(y ~ jit, data = df, col = col, pch = pch)
	return (df)
}

multibox = function (data, groups, ...) {
#	## directonality validation
#	young_age = groups$young_age
#	old_age = groups$old_age
#
#	age_groups = list(young = 0:(young_age), middle = (young_age+1):(old_age-1), old = old_age:69, midold=70:79, veryold=80:90)
#	labels = paste(sapply(age_groups, min), sapply(age_groups, max), sep='-')
#	par(mfrow=c(4, 4))
#	m.tum.logit = log2(as.matrix(m.tum)) - log2(1 - as.matrix(m.tum))
#	p.tum = intersect(names(tu.ages(st)), colnames(m.tum.logit))
#	m.tum.logit = m.tum.logit[, p.tum]
#	a.tum = tu.ages(st)[p.tum]
#	age_group.tum = sapply( round(unname(a.tum)), function (x) which(sapply (sapply( age_groups , is.element, x), any)))
#	m.norm.logit = log2(as.matrix(m.norm)) - log2(1 - as.matrix(m.norm))
#	p.norm = intersect(names(tu.ages(st)), colnames(m.norm.logit))
#	m.norm.logit = m.norm.logit[, p.norm]
#	a.norm = tu.ages(st)[p.norm]
#	age_group.norm = sapply( round(unname(a.norm)), function (x) which(sapply (sapply( age_groups , is.element, x), any)))
#	i=1
#	if (de.meth.tumors[i, 'genenormal'] == T)
#		next
#	cpg = rownames(de.meth.tumors)[i]
#	gene = st$meth450genes[cpg, 'gene']
##				plot(as.numeric(m.logit[cpg, ]) ~ tu.ages(st)[colnames(m)], main = sprintf('lfc and methylation validation, %s, lfc=%.2f, qval=%.2e', cpg, de.meth.tumors$intercept[i], de.meth.tumors$qval[i]))
#	x.norm = age_group.norm + runif(length(age_group.norm))/2 - 0.25 
#	x.tum  = age_group.tum  + runif(length(age_group.tum))/2 - 0.25 
#	boxplot(as.numeric(m.tum.logit[cpg, ]) ~ age_group.tum, names=labels, main = sprintf('meth,%s,%s, int=%.3f, q=%.2f, nl.q=%.3f ', cpg, gene, de.meth.tumors$intercept[i], de.meth.tumors$qval[i], de.meth.normals[cpg, 'qval']), cex.main = 0.5, cex.axis=0.5, outline=F)
#	points(as.numeric(m.norm.logit[cpg, ]) ~ x.norm, col=rgb(0.0,0.0,0.4,0.3), pch=20)
#	points(as.numeric(m.tum.logit[cpg, ]) ~ x.tum, col=rgb(1,0,0,0.2), pch=20)
} 


