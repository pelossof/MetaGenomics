source('~/paths.R')
source.util('de')


read.signatures = function (sigpath = '/Users/pelossof/Projects/colorectal/kras/signatures') {
	dtf = read.csv(paste(sigpath,'dtf.csv', sep='/'), stringsAsFactors=F)
	hacohen = read.csv(paste(sigpath,'hacohen.csv', sep='/'), stringsAsFactors=F, header=F)
	immgen.macrophage = read.csv(paste(sigpath,'immgen.macrophage.csv', sep='/'), header=T, stringsAsFactors=F);
	csf1r = read.csv(paste(sigpath,'csf1r.csv',sep='/'), stringsAsFactors=F, header=F)
	estimate = read.csv(paste(sigpath,'estimate-genes.csv',sep='/'), stringsAsFactors=F, header=F, row.names=1)
	est.imm = unlist(toupper(estimate['Immune141_UP', 2:ncol(estimate)]))
	est.str = unlist(toupper(estimate['Stromal141_UP', 2:ncol(estimate)]))
	load(paste(sigpath,'Galon_Engler_list.RData', sep='/'))

	dtf = toupper(unlist(dtf))
	hacohen.all = toupper(unlist(hacohen$V2))
	mph.immgen = toupper(unlist(immgen.macrophage$Gene.Symbol))
	csf1r = toupper(unlist(csf1r))

	signatures.list = list (stroma = list(dtf = dtf, est.str = est.str), macrophage = list(est.imm = est.imm, hacohen = hacohen.all, mph = mph.immgen, csf1r = csf1r), 
		Tcells = c( galon_engler[grep ('CD8', names(galon_engler), value=T)], galon_engler[grep ('CD4', names(galon_engler), value=T)]),  
		galon_engler[sapply(galon_engler, length)>=50])
	return (signatures.list)
}

plot.signatures = function(de.table, filename, signatures.list = read.signatures(),  lfc = "log2FoldChange", padj = "padj", sig.fdr = 0.2, ...) {
	pdf(filename)#, width=10, height=5)
	for (signatures in signatures.list) {
		for (s in 1:length(signatures)) {
			volcano.plot(de.table, signatures[[s]], lfc.name = lfc, padj.name = padj, main=sprintf('%s', names(signatures)[s]), arrows=T,...)
		}
	}
	dev.off()

}
