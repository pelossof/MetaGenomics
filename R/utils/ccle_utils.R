source('~/paths.R')


library(preprocessCore)
library('hgu133plus2.db')
source(getUtilsFname('affy_utils.R'))


fCCLE.loadCOAD <- function (path = "CCLE_Expression.Arrays_2013-03-18", ann.file = "CCLE_sample_info_file_2012-10-18.txt", skipLoad = FALSE) {
	if (skipLoad) {
		ccle.ann <- read.table ( "CCLE_sample_info_file_2012-10-18.txt" ,  header = T, sep = '\t', stringsAsFactors = F)
		f.nofile <- which(ccle.ann$Expression.arrays=='')
		ccle.ann <- ccle.ann[-f.nofile, ]
		coad.ix <- grep ( 'LARGE_INTESTINE', ccle.ann$CCLE.name )
#		coad.ix <- coad.ix[1:3]
		coad.names <- fCCLE.getname(ccle.ann$Cell.line.primary.name[coad.ix])
		f <- ccle.ann$Expression.arrays[coad.ix]
		coad.files <- paste(path, paste(f,'.CEL',sep=''), sep = '/')
		cels <- read.cel.files(coad.files)
		cels.norm <- normalize.cel.data(cels)
		colnames(cels.norm) <- coad.names

		save(cels.norm, file = getDataFname('processed/ccle.loadCOAD.rd'));
	} else {
		load(getDataFname('processed/ccle.loadCOAD.rd'));
	}
	return (cels.norm)
}

fCCLE.readmaf <- function (file="CCLE_Ann/CCLE_hybrid_capture1650_hg19_NoCommonSNPs_NoNeutralVariants_CDS_2012.05.07.maf") {
	file <- getDataFname(file)
	coad.maf <- read.delim(file, header = T, sep = '\t', stringsAsFactors = F)
	maf <- cbind (Hugo_Symbol=coad.maf[,'Hugo_Symbol'], Tumor_Sample_Barcode=fCCLE.getsample(coad.maf[,'Tumor_Sample_Barcode']))
	return (maf)
}

fCCLE.getsample <- function (barcodes) {
	barcodes.list <- strsplit (barcodes , split = '_')
	return (fCCLE.getname(sapply( barcodes.list, '[[',1) ))
}

fCCLE.getname <- function (cel.names) {
	return (toupper(gsub('[-, ,\\.]','',cel.names)))
}


fCCLE.getmuts <- function (maf, genes) {
	#get kras, braf, tp53, apc mutations

	muts <- matrix(0, 4, length( unique( maf [,2])))
	colnames(muts) <- unique(maf[,2])
	rownames(muts) <- c('KRAS', 'BRAF','TP53', 'APC')

	kras <- maf[grep ('KRAS', maf[,1]),2]
	braf <- maf[grep ('BRAF', maf[,1]),2]
	tp53 <- maf[grep ('TP53', maf[,1]),2]
	apc <- maf[grep ('APC', maf[,1]),2]

	muts['KRAS',kras] <- 1
	muts['BRAF',braf] <- 1
	muts['TP53',tp53] <- 1
	muts['APC' ,apc] <- 1

	return(t(muts))
}

