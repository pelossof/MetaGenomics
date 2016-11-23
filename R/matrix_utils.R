## matrix_utils.R


reindex = function(tab, index) {
# rename rownames of a matrix with a new index
# rows that are duplicates of the index or where the index
# is NA re removed
# then the index becomes the new rownames of tab
	ix = c()
	if (length(index) == 1)
		if (index %in% colnames(tab)) ix = tab[,index]
	if (length(index) == nrow(tab)) ix = index
	if(length(ix) == 0) {
		print('Error in reindex, returning original table')
		return (tab)
	}

	ix.unique = !duplicated(ix) & !is.na(ix)
	tab = tab[ix.unique, ]
	rownames(tab) = ix[ix.unique]
	return (tab)
}

intersect.mats = function(x, y, byRow = T, byCol = T) {
# takes two matrices x,y and intersects their col and row names
# returns a list with new matrices with intersecting indexes
	col.ix.x = intersect(colnames(x), colnames(y))
	row.ix.x = intersect(rownames(x), rownames(y))
	col.ix.y = col.ix.x
	row.ix.y = row.ix.x

	if (byRow == F) {
		row.ix.x = rownames(x)
		row.ix.y = rownames(y)
	}

	if (byCol == F) { 
		col.ix.x = colnames(x)
		col.ix.y = colnames(y)	
	}
	
	return (list(x[row.ix.x, col.ix.x], y[row.ix.y, col.ix.y]))
}



jga2tcga_barcode = function (ids) {
	project = 'JGA'
	tss = 'DC' 		# MSKCC Rectum adenocarcinoma
	ids = gsub('_|-|\\.','',ids)
	ids = toupper(ids)
	ids_parse = str_match(ids, '(JM|AS|TM|TS)([0-9]+)([BTN]+)([0-9]*)')
	participant = paste0(ids_parse[,2], ids_parse[,3])
	btn = ids_parse[,4]
	sample = rep(NA, length(btn))
	sample[grep('B', btn)] = '01' #01 primary solid tumor
	sample[grep('N', btn)] = '11' #11 solid tissue normal
	sample[grep('T', btn)] = '02' #02 reccurent tumor
	portion = rep("1", nrow(ids_parse))
	portion[ids_parse[,5] != ""] = ids_parse[ids_parse[,5] != "",5]
	return (paste(project, tss, participant, sample, portion, sep = '-'))
}

res2gene = function (res) {
	res$gene = sapply(strsplit(rownames(res), '|', fixed = T), '[[', 1)
	res = subset(res, !is.na(res$gene) & !duplicated(res$gene))
	rownames(res) = res$gene
	res = res[, -which(colnames(res) == 'gene')]
	return (res)
}

siglist = function(res, study = "", padj = 0.05, padj.name = 'padj', lfc.name = 'log2FoldChange', base.name = 'baseMean', base.min = 0, lfc.abs.min = 0) {
	res_clean = res[!is.na(res[,padj.name]) & !is.na(res[,lfc.name]) & !is.na(res[,base.name]) & !is.na(rownames(res)),   ]
	sig = res_clean[,padj.name]
	lfc = res_clean[, lfc.name]
	basemean = res_clean[, base.name]
	up = rownames(res_clean)[(sig <= padj) & (lfc >=  lfc.abs.min) & (basemean >= base.min) & !is.na(padj) & !is.na(rownames(res_clean))]
	dn = rownames(res_clean)[(sig <= padj) & (lfc <= -lfc.abs.min) & (basemean >= base.min) & !is.na(padj) & !is.na(rownames(res_clean))]
	ret = list(up, dn)
	if (study != "") study = paste0(study, '_')
	names(ret) = c(paste0(study, 'UP'), paste0(study, 'DOWN')) 
	return(ret)
}


