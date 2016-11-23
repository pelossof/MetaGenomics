maf2mat <- function (maf, gene.col = 'Hugo_Symbol', id.col = 'Tumor_Sample_Barcode') {
	maf <- cbind (maf[,gene.col], maf[,id.col])

	genes <- unique(maf[ , 1]) 
	pats = unique(maf[ , 2])
	mafmat = matrix(0, nrow = length(genes), ncol = length(pats))
	rownames(mafmat) <- genes
	colnames(mafmat) <- pats
	for (p in pats) {
		g <- unique(maf[maf[,2] %in% p, 1])
		mafmat[g, p] <- 1
	}

	return (mafmat)
}


impute.mat = function (mat, impute = list(genes = NA, val = NA)) {
# example
#B = matrix( c(2, 4, 3, 1, 5, 7), nrow=3, ncol=2, dimnames=list(c('a', 'b','c'), c('1', '2')))
#impute.mat(B, list(list(val = 0, genes = c('d', 'e', 'f', 'c')), list(val = NA, genes = c('d', 'e', 'g', 'h'))))  
#
	mat.new = mat
	for (imp in impute) {
		genes = imp$genes
		val = imp$val
		genes.original = rownames(mat.new)
		genes.impute = setdiff(genes, rownames(mat.new))
		mat.new = rbind(mat.new, matrix(val, length(genes.impute), ncol(mat.new)))
		rownames(mat.new) = c(genes.original, genes.impute)
	}
	return (mat.new)
}

load.impact = function () {
	impact = list()
	impact[[1]] = toupper(unname(unlist(read.csv(paths.data('impact/impact230.csv'), stringsAsFactors = F, header = F))))
	impact[[2]] = toupper(unname(unlist(read.csv(paths.data('impact/impact341.csv'), stringsAsFactors = F, header = F))))
	impact[[3]] = toupper(unname(unlist(read.csv(paths.data('impact/impact410.csv'), stringsAsFactors = F, header = F))))
	return(impact)
}

