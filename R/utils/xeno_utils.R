require(org.Hs.eg.db)
require(org.Mm.eg.db)

source('~/paths.R')

mm2hg.ens.map = read.csv(paths.data('xeno/mm10ENS_hg19GENE.csv'), header= T, stringsAsFactors = F)
mm2hg.ens.map = mm2hg.ens.map[!duplicated(mm2hg.ens.map[,1]),]
rownames(mm2hg.ens.map) = mm2hg.ens.map[, 1]
colnames(mm2hg.ens.map) = c('ensid', 'mm', 'hg')


mm2hg.ens = function (ensids) {
	hg = mm2hg.ens.map[ensids,'hg']
	names(hg) = ensids
	return(hg)
}

ref2ez <- function (refid, db = 'Hs') {
	# refseq to entrez
	if (db == 'Hs')
		x <- org.Hs.egREFSEQ2EG
	else
		x <- org.Mm.egREFSEQ2EG

	mapped_seqs <- mappedkeys(x)
	xx <- as.list(x[mapped_seqs])

	ezid = unlist(xx[refid])
	return (ezid)	
}

ez2sym <- function (ezid, db = 'Hs') {
	# entrez to symbol 
	if (db == 'Hs')
		x <- org.Hs.egSYMBOL
	else
		x <- org.Mm.egSYMBOL
	mapped_genes <- mappedkeys(x)
	xx <- as.list(x[mapped_genes])

	symbol <- unlist(xx[ezid])
	return (symbol)
}


