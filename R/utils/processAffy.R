# processSamples.R
# Raphael Pelossof 2013

## process all samples
# copy number through analyzeCGH
# gene expression through analyzeAffy
# microRNA through analyzeMir
# point mutations through analyzeMut

source ('~/paths.R' )

source(getUtilsFname("affy_utils.R"))
library("hgu133plus2.db")
require(glmnet)
require(limma)
require(gplots)
require(preprocessCore)
require(stringr)
#libraries for go analysis
library("org.Hs.eg.db")


loadTiming <- function ( skipLoad = FALSE, mas5.thresh = 0.85 ) {
	if(skipLoad) {
		file.remove(getProcessedFname('timing.norm.Rdata'))
		file.remove(getProcessedFname('timing.biopsy.Rdata'))
		file.remove(getProcessedFname('timing.tumor.Rdata'))

		#load normals
		filenames <- list.files(getDataFname('timing.cel'), pattern = "*.CEL", full.names = TRUE)
		filenames <- filenames[!grepl ('AS', filenames)]  #remove aoosog samples
		cel.names <- str_match(filenames, '([0-9]+)N')[,2]
		cel.names <- paste(cel.names,'-N',sep="")
		if(!file.exists(getProcessedFname('timing.norm.Rdata'))) {
			cels.n <- processCels(filenames, mas5.thresh, cel.names)
			save(cels.n, file=getProcessedFname('timing.norm.Rdata'))
		}

		#load biopsies
		filenames <- list.files("data/mRNA/R01_Biopsy_Before_Treatment", pattern = "*.CEL", full.names = TRUE)
		filenames <- filenames[!grepl ('AS', filenames)]  #remove aoosog samples
		cel.names <- str_match(filenames, '([0-9]+)B')[,2]
		cel.names <- paste(cel.names,'-B',sep="")
		if(!file.exists("data/process/timing.biopsy.Rdata")) {
			cels.b <- processCels(filenames, mas5.thresh, cel.names)
			save(cels.b, file="data/process/timing.biopsy.Rdata")
		}

		#load tumors
		filenames <- list.files("data/mRNA/R01_Tumor_After_Treatment", pattern = "*.CEL", full.names = TRUE)
		filenames <- filenames[!grepl ('AS', filenames)]  #remove aoosog samples
		cel.names <- str_match(filenames, '([0-9]+)T')[,2]
		cel.names <- paste(cel.names,'-T',sep="")
		if(!file.exists("data/process/timing.tumor.Rdata")) {
			cels.t <- processCels(filenames, mas5.thresh, cel.names)
			save(cels.t, file="data/process/timing.tumor.Rdata")
		}

	} else {
		load(getProcessedFname('timing.norm.Rdata'))
		load(getProcessedFname('timing.biopsy.Rdata'))
		load(getProcessedFname('timing.tumor.Rdata'))
	}


	return (list (clin = loadTimingClin(), mut = loadTimingMuts(), 
		arr = cbind(cels.n, cels.b, cels.t), arr.type="hgu133plus2.db"))
}

loadTimingClin <- function () {
	clin <- read.table(getDataFname("timing/ListOfAllSamples1.csv"),header=TRUE, row.names="ID", sep=",",stringsAsFactors=F)
	ix <- !grepl ('AS', rownames(clin))
	clin <- clin[ix, ]

	#name a few columns
	colnames(clin)[which (colnames(clin)=="pCR.N.pCR")] <- 'pcr'
	colnames(clin)[which (colnames(clin)=="TRG")] <- 'trg'

	clin[clin=='NA'] <- NA
	clin[clin=='?'] <- NA
	clin[clin=='Y'] <- 1
	clin[clin=='N'] <- 0
	#clin[clin=='I'] <- 1
	#clin[clin=='II'] <- 2
	#clin[clin=='III'] <- 3

	return (clin)
}

loadTimingMuts <- function () {
	data <- read.table(getDataFname("timing/ListOfAllSamples1.csv"),header=TRUE, row.names="ID", sep=",", stringsAsFactors=F)
	ix <- !grepl ('AS', rownames(data))
	data <- data[ix, ]

	data[data == 'NA'] <- NA
	data[data == '?'] <- NA

	n <- nrow(data) #ncol(data.probes)
	colnames <- c('krasmut','kraswt','kras12','kras13','tp53','braf')
	info <- matrix (rep (0, length(colnames)*n), ncol=length(colnames), 
		dimnames=list(rownames(data), colnames))
	info[,"krasmut"] <- data[,"Kras"]!="wt"
	info[,"kraswt"] <- data[,"Kras"]=="wt"
	info[,"kras12"] <- grepl("12", data[,"Kras"])
	info[,"kras13"] <- grepl("13", data[,"Kras"])
	info[,"tp53"] <- data[,"P53"]!="wt"
	info[,'braf'] <- !(data[,'Braf'] =='wt' | data[,'Braf']=='?')
#	info[,'PIK3CA'] <- data[,'PIK3CA']!='wt'
	return (info)
}

loadAcosog <- function ( skipLoad = FALSE, mas5.thresh = 0.85 ) {
	#load biopsies
	filenames <- list.files("data/mRNA/R01_Biopsy_Before_Treatment", pattern = "*.CEL", full.names = TRUE)
	filenames <- filenames[grepl ('AS', filenames)]  #remove aoosog samples
	cel.names <- str_match(filenames, '(AS[0-9]+)B')[,2]
	cel.names <- paste(cel.names,'-B',sep="")
	if(skipLoad) {
		file.remove ("data/process/acosog.biopsy.Rdata")
	} else {
		load ("data/process/acosog.biopsy.Rdata")
	}

	if(!file.exists("data/process/acosog.biopsy.Rdata")) {
		cels.b <- processCels(filenames, mas5.thresh, cel.names)
		save(cels.b, file="data/process/acosog.biopsy.Rdata")
	}

	return (list (clin = loadAcosogClin(), mut = loadAcosogMuts(), 
		arr = cels.b, arr.type="hgu133plus2.db"))
}

loadAcosogClin <- function () {
	clin <- read.table("data/ListOfAllSamples1.csv",header=TRUE, row.names="ID", sep=",",stringsAsFactors=F)
	ix <- grepl ('AS', rownames(clin))
	clin <- clin[ix, ]

	#name a few columns
	colnames(clin)[which (colnames(clin)=="pCR.N.pCR")] <- 'pcr'
	colnames(clin)[which (colnames(clin)=="TRG")] <- 'trg'

	clin[clin=='NA'] <- NA
	clin[clin=='?'] <- NA
	clin[clin=='Y'] <- 1
	clin[clin=='N'] <- 0
	#clin[clin=='I'] <- 1

	return (clin)
}

loadAcosogMuts <- function () {
	data <- read.table("data/ListOfAllSamples1.csv",header=TRUE, row.names="ID", sep=",")
	ix <- grepl ('AS', rownames(data))
	data <- data[ix, ]
	n <- nrow(data) #ncol(data.probes)
	colnames <- c('krasmut','kraswt','kras12','kras13','tp53','braf')
	info <- matrix (rep (0, length(colnames)*n), ncol=length(colnames), 
		dimnames=list(rownames(data), colnames))
	info[,"krasmut"] <- data[,"Kras"]!="wt"
	info[,"kraswt"] <- data[,"Kras"]=="wt"
	info[,"kras12"] <- grepl("12", data[,"Kras"])
	info[,"kras13"] <- grepl("13", data[,"Kras"])
	info[,"tp53"] <- data[,"P53"]!="wt"
	info[,'braf'] <- !(data[,'Braf'] =='wt' | data[,'Braf']=='?')
#	info[,'PIK3CA'] <- data[,'PIK3CA']!='wt'
	return (info)
}


processCels <- function (filenames, mas5.thresh = 0.85, names = NA) {
#set mas5.thresh=NA to skip mas5 probe removal	
	cels <- read.cel.files(filenames)

	#normalize probes with expresso
	norm.probes <- normalize.cel.data(cels)
	if (missing(names)) {
		colnames(norm.probes) <- filenames
	} else {
		colnames(norm.probes) <- names
	}

	#use mas5 to remove lowely lit probes
#	if (!is.na (mas5.thresh)) {
#		m5.calls <- mas5calls(cels)
#		m5.calls.x <- exprs(m5.calls)
#		absent.probes <- apply(m5.calls.x=="A", 1,sum) > mas5.thresh * ncol(m5.calls.x)
#		norm.probes <- norm.probes[!absent.probes, ]
#	}		
	return (norm.probes)
}


loadTCGA <- function (skipLoad = FALSE) {
#load tcga level 3 data
	tcga.filenames <- list.files(path="data/mRNA/tcga/Expression-Genes/UNC__AgilentG4502A_07_3/Level_3/", pattern="*", full.names = TRUE)
	tcga.ids <- substr(tcga.filenames, 95, 106)

	x <- lapply(tcga.filenames, function(f) {
		y <- readLines(f)[-1]
		as.numeric(sapply(strsplit(y, '\t'), '[', 3))
	})
	tcga.norm <- as.matrix(do.call(cbind, x))
	tcga.norm <- normalize.quantiles (tcga.norm)
	
	colnames(tcga.norm) <- tcga.ids
	y <- readLines(tcga.filenames[1],-1)
	rownames(tcga.norm) <- sapply(strsplit(y[2:length(y)], '\t'), '[', 2)



	#note, no intersection between sample names needed since validate by hand
	tcga.kras <- read.table("data/mRNA/tcga/mut/kras.txt", row.names = 1, header=TRUE,sep="\t")
	tcga.tp53 <- read.table("data/mRNA/tcga/mut/tp53.txt", row.names = 1, header=TRUE,sep="\t")
	n <- length(unique(tcga.ids))
	tcga <- matrix (rep (0, 5*n), ncol=5, dimnames=list(unique (tcga.ids), c('kraswt','krasmut','kras12','kras13','tp53')))
	tcga[rownames(tcga.kras),"krasmut"] <- 1
	tcga[setdiff(colnames(tcga.norm), rownames(tcga.kras)),"kraswt"] <- 1
	tcga[rownames(tcga.tp53),"tp53"] <- 1
	tcga[rownames(tcga.kras)[grep("13", tcga.kras[,"Mut"])], "kras13"] <- 1
	tcga[rownames(tcga.kras)[grep("12", tcga.kras[,"Mut"])], "kras12"] <- 1

	#remove replicat samples, possibly better average dups 
	dups <- duplicated (tcga.ids)
	tcga.norm <- tcga.norm[, !dups]
#	tcga.ag <- aggregate (t(tcga.norm), by=list (tcga.ids), FUN=mean)
#	tcga.norm1 <- t(tcga.ag)
	return(list(info = tcga, arr = tcga.norm))
}



loadTCGA2 <- function (skipLoad = FALSE) {
#load tcga level 2 data
	tcga.path <- "data/mRNA/tcga-level2"
	tcga.manifest <- read.table(paste(tcga.path, 'file_manifest.txt',sep='/'), 
		sep='\t', header = TRUE, row.names = "File.Name")
	tcga.manifest <- tcga.manifest[-c(1,2), ]
	rownames(tcga.manifest) <- tcga.manifest[,'File.Name']
	
	tcga.filenames <- list.files(path=paste(tcga.path, 'Expression-Genes/UNC__AgilentG4502A_07_3/Level_2' ,sep='/'), pattern="*", full.names = TRUE)
#	tcga.ids <- substr(tcga.filenames, 127, 138)

	x <- lapply(tcga.filenames, function(f) {
		y <- readLines(f)[-1]
		as.numeric(sapply(strsplit(y, '\t'), '[', 2))
	})
	tcga.norm <- as.matrix(do.call(cbind, x))
	tcga.norm <- tcga.norm[-1,]
	tcga.norm <- normalize.quantiles (tcga.norm)
	
	files <- sapply(strsplit(tcga.filenames, '/'), '[', 12) #get the filename
	tcga.ids <- as.character(tcga.manifest[files, 'Sample'])
	colnames(tcga.norm) <- tcga.ids

	y <- readLines(tcga.filenames[1],-1)
	rownames(tcga.norm) <- sapply(strsplit(y[3:length(y)], '\t'), '[', 1)



	#note, no intersection between sample names needed since validate by hand
	tcga.kras <- read.table("data/mRNA/tcga/mut/kras.txt", row.names = 1, header=TRUE,sep="\t")
	tcga.tp53 <- read.table("data/mRNA/tcga/mut/tp53.txt", row.names = 1, header=TRUE,sep="\t")
	n <- length(unique(tcga.ids))
	tcga <- matrix (rep (0, 5*n), ncol=5, dimnames=list(unique (tcga.ids), c('kraswt','krasmut','kras12','kras13','tp53')))
	tcga[rownames(tcga.kras),"krasmut"] <- 1
	tcga[setdiff(colnames(tcga.norm), rownames(tcga.kras)),"kraswt"] <- 1
	tcga[rownames(tcga.tp53),"tp53"] <- 1
	tcga[rownames(tcga.kras)[grep("13", tcga.kras[,"Mut"])], "kras13"] <- 1
	tcga[rownames(tcga.kras)[grep("12", tcga.kras[,"Mut"])], "kras12"] <- 1

	#remove replicat samples, possibly better average dups 
	dups <- duplicated (tcga.ids)
	tcga.norm <- tcga.norm[, !dups]
#	tcga.ag <- aggregate (t(tcga.norm), by=list (tcga.ids), FUN=mean)
#	tcga.norm1 <- t(tcga.ag)
	return(list(info = tcga, arr = tcga.norm))
}

loadTCGA3 <- function (skipLoad = FALSE) {
	if(skipLoad) {
		file.remove("data/process/tcga3.Rdata")
	} else {
		load("data/process/tcga3.Rdata")
	}

	tcga.path <- "data/mRNA/tcga/colorectal-expr"
	arr.path <- "/Expression-Genes/UNC__AgilentG4502A_07_3/Level_3"
	tcga.manifest <- read.table(paste(tcga.path, 'file_manifest.txt',sep='/'), 
		sep='\t', header = TRUE, row.names = "File.Name")
	tcga.manifest <- tcga.manifest[-c(1,2), ]
	rownames(tcga.manifest) <- tcga.manifest[,'File.Name']
	
	tcga.filenames <- list.files(path=paste(tcga.path, 'Expression-Genes/UNC__AgilentG4502A_07_3/Level_2' ,sep='/'), pattern="*", full.names = TRUE)



	#load normals
	filenames <- list.files("data/mRNA/R01_Normal", pattern = "*.CEL", full.names = TRUE)
	filenames <- filenames[!grepl ('AS', filenames)]  #remove aoosog samples
	cel.names <- str_match(filenames, '([0-9]+)N')[,2]
	cel.names <- paste(cel.names,'-N',sep="")
	if(!file.exists("data/process/timing.norm.Rdata")) {
		cels.n <- processCels(filenames, mas5.thresh, cel.names)
		save(cels.n, file="data/process/timing.norm.Rdata")
	}

	#load biopsies
	filenames <- list.files("data/mRNA/R01_Biopsy_Before_Treatment", pattern = "*.CEL", full.names = TRUE)
	filenames <- filenames[!grepl ('AS', filenames)]  #remove aoosog samples
	cel.names <- str_match(filenames, '([0-9]+)B')[,2]
	cel.names <- paste(cel.names,'-B',sep="")
	if(!file.exists("data/process/timing.biopsy.Rdata")) {
		cels.b <- processCels(filenames, mas5.thresh, cel.names)
		save(cels.b, file="data/process/timing.biopsy.Rdata")
	}

	#load tumors
	filenames <- list.files("data/mRNA/R01_Tumor_After_Treatment", pattern = "*.CEL", full.names = TRUE)
	filenames <- filenames[!grepl ('AS', filenames)]  #remove aoosog samples
	cel.names <- str_match(filenames, '([0-9]+)T')[,2]
	cel.names <- paste(cel.names,'-T',sep="")
	if(!file.exists("data/process/timing.tumor.Rdata")) {
		cels.t <- processCels(filenames, mas5.thresh, cel.names)
		save(cels.t, file="data/process/timing.tumor.Rdata")
	}

	return (list (clin = loadTimingClin(), mut = loadTimingMuts(), 
		arr = cbind(cels.n, cels.b, cels.t), arr.type="hgu133plus2.db"))
}




loadReid <- function (skipLoad=FALSE) {
#load reid data
	files <- list.files(path = "data/mRNA/E-GEOD-20842.processed.1", pattern = "GSM*", full.names = TRUE)
	ag<- lapply(files, read.table, sep="\t", header=TRUE, row.names="Reporter.Identifier")
	reid.norm <- as.matrix(do.call(cbind, ag))
	reid.ids <- substr(files,36,44)
	colnames(reid.norm) <- reid.ids

	darkprobe <- min(reid.norm["DarkCorner",])
	reid.max <- apply(reid.norm, 1, max)
	reid.norm <- reid.norm[reid.max > darkprobe,] #keep only probes where the darkest probe is bright enough
	reid.genes <- unlist(mget(rownames(reid.norm),hgug4112aSYMBOL, ifnotfound=NA))
	reid.norm <- reid.norm[!is.na(rownames(reid.norm)),]

	reid.an <- read.table("data/mRNA/E-GEOD-20842.processed.1/E-GEOD-20842.sdrf.txt", header=TRUE,sep="\t", row.names="Hybridization.Name")
	stopifnot(length(rownames(reid.an)) == length(unique(rownames(reid.an))))
	cols <- c('id','male','krasmut','kraswt','mucosa','tumor')
	n <- length(rownames(reid.an))
	m <- length(cols)
	reid <- matrix (rep (0, m*n), ncol=m, dimnames=list(rownames(reid.an), cols))
	reid[,'id'] <- as.numeric(substr(as.character(reid.an[,"FactorValue..PATIENT.ID."]),2,5))
	reid[,"male"] <- reid.an[,"FactorValue..GENDER."]=="male"
	reid[,"krasmut"] <- grepl("mut",reid.an[,"FactorValue..GENOME.VARIATION."])
	reid[,"kraswt"] <- grepl("wild",reid.an[,"FactorValue..GENOME.VARIATION."])
	reid[,"mucosa"] <- reid.an[,"FactorValue..TISSUE."]=="mucosa"
	reid[,"tumor"] <- !reid.an[,"FactorValue..TISSUE."]=="mucosa"

	return(list(info = reid, arr = reid.norm))
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

	if (nrow(tt)>0)
		tt$gene <- probes.to.symbols(tt$ID, arr.type)

	return (tt)
}


aggregateArray <- function (ar, chip, collapse = FALSE) {
#ar is a probe intensity matrix where each row is a probe, and each column is a sample
	if (chip == "hgug4112aSYMBOL") { #agilent
		ar.genes <- unlist(mget(rownames(ar),hgug4112aSYMBOL, ifnotfound=NA))
	}
	else {
		ar.genes <- probes.to.symbols(rownames(ar), "hgu133plus2.db")
	}	
	ar.ag <- aggregate (ar, by=list (ar.genes), FUN=median)
	rownames(ar.ag) <- ar.ag[[1]]
	ar.ag <- ar.ag[,-1]
	if (collapse)
		ar.ag <- apply(ar.ag,1,median)

	return (ar.ag)
}
