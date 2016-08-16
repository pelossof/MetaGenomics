library(parallel)
library(preprocessCore)
library(RCurl)
library(rjson)
source('~/paths.R')
source.util('survival')

tu.uuid2barcode <- function(uuids) {
	# Query TCGA's UUID to barcode Web Service.
	resp <- getURL("https://tcga-data.nci.nih.gov/uuid/uuidws/mapping/json/uuid/batch", customrequest="POST", httpheader=c("Content-Type: text/plain"), postfields=paste(uuids, collapse=","))
	 
	# Extract mappings from response.
	mappings <- fromJSON(resp)$uuidMapping
	m <- unlist(mappings)
	uu2bar <- m[seq(1,length(m), 2)]
	names(uu2bar) = m[seq(2,length(m), 2)]
	return (uu2bar)
}

parse.name.prot <- function (filename) {
	strbits <- strsplit(basename(filename),'.', fixed=T)
	return (strbits[[1]][6])
}

parse.name.expr <- function (filename) {
	strbits <- strsplit(basename(filename),'.', fixed=T)
	return (strbits[[1]][3])
}

parse.name.tcga <- function (filename) {
	fname <- basename(filename)
	ix <- regexpr('TCGA', fname)[1]
	return(substr(fname, ix, ix+27))
}

tu.tumors <- function (mat, simplify.barcode = FALSE) {
	mat <- mat[, -c(which(!tu.istumor(colnames(mat))), ncol(mat)+1)]
	if (simplify.barcode) {
		colnames(mat) <- tu.getsample(colnames(mat))
	}
	return(mat)
}

tu.normals <- function (mat, simplify.barcode = FALSE) {
	mat <- mat[, -c(which(!tu.isnormal(colnames(mat))), ncol(mat)+1)]
	if (simplify.barcode) {
		colnames(mat) <- tu.getsample(colnames(mat))
	}
	return(mat)
}

tu.ages <- function(study) {
	ages = study$clin$patient$days_to_birth/-365
	names(ages) = study$clin$patient$bcr_patient_barcode
	ages = ages[!is.na(ages)]
	return (ages)
}


tu.genes <- function (mat) {
	return (unlist(lapply(strsplit(rownames(mat), '|', fixed = T),'[[',1)))
}

tu.ctsymbols = function (mat) { 
	gene.names = tu.genes(mat)
	ix = !duplicated(gene.names) & !is.na(gene.names)
	mat = mat[ix, ]
	rownames(mat) = gene.names[ix]
	return (mat)
}

tu.plotkm = function(st, groups, ...) {
	death = !is.na(st$clin$patient$days_to_death)
	time = st$clin$patient$days_to_last_followup
	names(time) = st$clin$patient$bcr_patient_barcode
	names(death) = st$clin$patient$bcr_patient_barcode
	time[is.na(time)] = st$clin$patient$days_to_death[is.na(time)]
	survu.plotkm(time, death, groups)
}


tu.lvl3 <- function (path, pattern = "*", gene.column = 2, numeric.column = 3, name.parse.func = parse.name.expr, skip=1, mc.cores=6) {
	# expecting files with: ID, Gene, Expr
	# if there is no ID use gene.column = 1, numeric.column = 2
	# expecting file name to be in the format: path/TCGA-ID.extension
	
	tcga.filenames <- list.files(path, pattern=pattern, full.names = TRUE)
	if(length(tcga.filenames) == 0) warning('No files found');
	#tcga.names <- sapply(strsplit(tcga.filenames, split = "__"), '[[', 4)
	#tcga.names <- sub("^([^.]*).*", "\\1", basename(tcga.filenames ))
	tcga.names <- unlist(lapply(tcga.filenames, name.parse.func))
	x <- mclapply(tcga.filenames, function(f) {
			y <- readLines(f)
			as.numeric(sapply(strsplit(y, '\t'), '[', numeric.column))
		}, mc.cores=mc.cores)		
	tcga.norm <- as.matrix(do.call(cbind, x))
	
	colnames(tcga.norm) <- tcga.names
	y <- readLines(tcga.filenames[1])
	rownames(tcga.norm) <- sapply(strsplit(y, '\t'), '[', gene.column)

	if (skip > 0)  {
		tcga.norm <- tcga.norm[-seq(1:skip), ]
	}

	return ( tcga.norm )
}


tu.getsample <- function (barcodes) {
	barcodes.list <- strsplit (barcodes , split = '-')
	return (sapply( barcodes.list, function (x) {paste(x[1:3],sep='-', collapse="-")})) 
}


tu.istumor <- function (barcodes) {
	barcodes.list <- strsplit (barcodes , split = '-')
	barcode.sample <- sapply( barcodes.list, '[', 4)
	barcode.number <- as.numeric(substr(barcode.sample, 1, 2))

	return (barcode.number < 10)
}


tu.isnormal <- function (barcodes) {
	barcodes.list <- strsplit (barcodes , split = '-')
	barcode.sample <- sapply( barcodes.list, '[', 4)
	barcode.number <- as.numeric(substr(barcode.sample, 1, 2))

	return (barcode.number > 10 & barcode.number < 20)
}

tu.issingle <- function (barcodes) {
	return ( grepl('single', barcodes))
}

tu.ispaired <- function (barcodes) {
	return ( grepl('paired', barcodes))
}

tu.readmaf <- function (file) {
	coad.maf <- read.delim(f, header = T, sep = '\t', stringsAsFactors = F)
	maf <- cbind (coad.maf[,'Hugo_Symbol'], tu.getsample(coad.maf[,'Tumor_Sample_Barcode']))
	return (maf)
}


tu.getmuts <- function (file, sample.names = NA) {
	#get kras, braf, tp53, apc mutations
	#sample names include all samples, also ones that may not appear in the maf file

	tcga.maf <- read.delim(file, header = T, sep = '\t', stringsAsFactors = F)
	maf <- cbind (tcga.maf[,'Hugo_Symbol'], tu.getsample(tcga.maf[,'Tumor_Sample_Barcode']))

	samples <- union( maf[,2] , tu.getsample(read.universe))
	muts <- matrix(0, 4, length( samples))
	colnames(muts) <- unique(samples)
	rownames(muts) <- c('KRAS', 'BRAF','TP53', 'APC')

	kras <- maf[grep ('KRAS', maf[,1]),2]
	braf <- maf[grep ('BRAF', maf[,1]),2]
	tp53 <- maf[grep ('TP53', maf[,1]),2]
	apc <- maf[grep ('APC', maf[,1]),2]

	muts['KRAS',kras] <- 1
	muts['BRAF',braf] <- 1
	muts['TP53',tp53] <- 1
	muts['APC' ,apc] <- 1

	return(muts)
}


tu.cut <- function (mat, row=NA, col=NA, simplify.barcode = F) {
	if (is.na(row)) row = rownames(mat)
	if (is.na(col)) col = colnames(mat)
	if (simplify.barcode) {
		colnames(mat) = tu.getsample(colnames(mat))
		names(col) = tu.getsample(names(col))
	}

	return (mat[row, col])
}


tu.getmutsall <- function (file, sample.names = NA) {
	#get kras, braf, tp53, apc mutations
	#sample names include all samples, also ones that may not appear in the maf file

	tcga.maf <- read.delim(file, header = T, sep = '\t', stringsAsFactors = F)
	maf <- cbind (tcga.maf[,'Hugo_Symbol'], tu.getsample(tcga.maf[,'Tumor_Sample_Barcode']))

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


tu.getmsi <- function (file) {
	# microsatellite-stable (MSS)
	# low level MSI (MSI-L) if less than 40% of markers were altered
	# high level MSI (MSI-H) if greater than 40% of markers were altered
	
	tcga.ms <- read.delim(f, header = T, sep = '\t', stringsAsFactors = F)
	ms <- tcga.ms[,'mononucleotide_and_dinucleotide_marker_panel_analysis_status']
	
	ms[ms=="MSS"] <- 0;
	ms[ms=="MSI-H"] <- 1;
	ms[ms=="MSI-L"] <- -1;
	ms <- as.integer(ms)
	names(ms) <- tcga.ms[,'bcr_patient_barcode']

	return(ms)
}

tu.read.counts <- function (study, library = 'pairedend') {
	mapfile <- paths.data(sprintf("TCGA/countslab.07_11_2014/samples_mapping_%s.txt", study))
	path <- paths.data(sprintf("TCGA/countslab.07_11_2014/tcga.%s.%s", study, library))
	print(path)
	print(mapfile)
	if(!file.exists(path)) {
		print(sprintf('could not find counts file %s/%s', study, library))
		return (NA)
	} else {
		print('found')
	}
	ct <- tu.lvl3(path, pattern = "*", gene.column = 1, numeric.column = 2, name.parse.func = function(x){basename(x)}, skip=0, mc.cores=6)
	ct.names <- gsub(sprintf('.cds.counts.%s.txt', library),'.bam', colnames(ct))
	map <- read.table (file = mapfile, stringsAsFactors = F, row.names = 3)
	colnames(map) <- c('uuid','barcode')
	if (is.na(map$barcode[1])) {
		# load barcodes from JSON
		map$barcode = tu.uuid2barcode(map$uuid)
	}

	colnames(ct) <- map[ct.names,'barcode']

	return (ct)
}

