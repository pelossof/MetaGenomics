source('paths.R')

## gdsc library
source('affy_utils.R')
library ('hgu133a.db')
source('ccle_utils.R')

fGDSC.loadCOAD <- function(path = '', skipLoad = FALSE) {
	ann.path = getDataFname('/gdsc')
	cels.path = getDataFname('/gdsc.cel')
	processed.path = getProcessedFname()
	
	if (skipLoad) { 
		#cell line name to array file name
		gdsc.files <- read.delim(paste(ann.path,'gdsc_files.txt',sep='/'), header = T, stringsAsFactors = F)

		#cell line name to cosmic id
		gdsc.lines <- read.csv(paste(ann.path,'gdsc_cell_lines_w2.csv',sep='/'), header = T, stringsAsFactors = F)

		#cosmic id to cell line type
		gdsc.ann <- read.delim(paste(ann.path,'CosmicCellLineProject_v65_280513.tsv',sep='/'), header=T, stringsAsFactors = F)

		gdsc.cos2site <- gdsc.ann$Primary.site
		names(gdsc.cos2site) <- gdsc.ann$ID_sample

		gdsc.cos2name <- gdsc.lines$CELL_LINE_NAME
		names(gdsc.cos2name) <- gdsc.lines$COSMIC_ID

		#get all the large intestine file names
		gdsc.coad.ix <- grep ('intestine', gdsc.cos2site)
		gdsc.coad.cos <- unique(names(gdsc.cos2site[gdsc.coad.ix]))

		gdsc.coad.names <- gdsc.cos2name[gdsc.coad.cos]
		gdsc.coad.files <- gdsc.files[gdsc.files$Source.Name %in% gdsc.coad.names,]

		fnames <- paste(cels.path, gdsc.coad.files$Array.Data.File,sep='/')

		cels <- read.cel.files(fnames)
		cels.norm <- normalize.cel.data(cels)
		colnames(cels.norm) <- fCCLE.getname(gdsc.coad.files$Source.Name)

		save(cels.norm, file = paste(processed.path, 'gdsc.loadCOAD.rd', sep = '/'));
	} else {
		load(file = paste(processed.path, 'gdsc.loadCOAD.rd', sep = '/'));
	}
	
	colnames(cels.norm) <- fCCLE.getname(colnames(cels.norm))
	return (list (mut=fGDSC.loadGDSCMuts(), arr = cels.norm, arr.type="hgu133a"))
}

fGDSC.loadGDSCMuts <- function(path = '', skipLoad = FALSE) {

	genes <- c('KRAS', 'BRAF', 'TP53', 'APC')
	gdsc.muts.all <- read.csv(getDataFname("gdsc/gdsc_mut.csv"),
			header=TRUE, row.names="Cell.Line", sep=",", stringsAsFactors=FALSE)
	gdsc.muts.genes <- gdsc.muts.all[,genes]
	gdsc.muts <- apply ( gdsc.muts.genes, c(1,2), function (x) grepl('p',x)) + 0
	rownames(gdsc.muts) <- fCCLE.getname(rownames(gdsc.muts))
	return (t(gdsc.muts))
}



