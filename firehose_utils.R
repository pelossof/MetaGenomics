#############################################
## Rafi Pelossof, MSK 2014
#
# firehose util - process pancan data for rbioportal
#


source('~/paths.R')
source(paths.utils('tcga_utils.R'))
source(paths.utils('maf_utils.R'))
source(paths.utils('cn_utils.R'))
source(paths.utils('clin_utils.R'))


errorFileNotFoundWarning <- function(e) {
	warning('File does not exist')
	return(NA)
}


fhu.version <- function(day, month, year) {
	return( data.frame(day=as.numeric(day), month=as.numeric(month), year=as.numeric(year)))
}

fhu.current.version = fhu.version(day = 18, month = 5, year = 2014)
fhu.current.version = fhu.version(day = 6, month = 12, year = 2014)
fhu.current.version = fhu.version(day = 1, month = 6, year = 2015)

analysis.current.version = fhu.version(day = 02, month = 04, year = 2015)


fhu.studies <- function (version = fhu.current.version) {
	day <- as.numeric(version$day)
	month <- as.numeric(version$month)
	year <- as.numeric(version$year)
	return(list.files(paths.data(sprintf("Firehose/stddata__%d_%.2d_%.2d", year, month, day))))
}

fhu.filename <- function(study, type, version = fhu.current.version) {
	day <- as.numeric(version$day)
	month <- as.numeric(version$month)
	year <- as.numeric(version$year)
	if (type == 'rppa') {
		type.str <- 'RPPA_AnnotateWithGene.Level_3'
		file <- sprintf('%s.rppa.txt', study)
	} else if (type == 'mut') {
		type.str <- 'Mutation_Packager_Calls.Level_3'
		file <- sprintf('%s.maf', study)
	} else if (type == 'mutraw') {
		type.str <- 'Mutation_Packager_Raw_Calls.Level_3'
		file <- sprintf('%s.maf', study)
	} else if (type =='clin') {
		#type.str <- 'Clinical_Pick_Tier1.Level_4'
#		file <- sprintf('%s.clin.stage2File.txt', study)
		type.str <- 'Merge_Clinical.Level_1'
		file <- sprintf('%s.clin.merged.txt', study)
	} else if (type == 'rna.ga') {
		type.str <- 'Merge_rnaseqv2__illuminaga_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.Level_3'
		file <- sprintf('%s.rnaseqv2__illuminaga_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.data.txt', study)
	} else if (type == 'rna.hi') {
		type.str <- 'Merge_rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.Level_3'
		file <- sprintf('%s.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.data.txt', study)
	} else if (type == 'rna.counts') {
		type.str <- 'Merge_rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes__data.Level_3'
		file <- sprintf('%s.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes__data.data.txt.counts.txt', study)
	} else if (type == 'cna') {
		type.str <- 'Merge_cna__illuminahiseq_dnaseqc__hms_harvard_edu__Level_3__segmentation__seg.Level_3'
		file <- sprintf('%s.cna__illuminahiseq_dnaseqc__hms_harvard_edu__Level_3__segmentation__seg.seg.txt', study)
	} else if (type == 'cnag') {
		type.str <- 'Merge_snp__genome_wide_snp_6__broad_mit_edu__Level_3__segmented_scna_minus_germline_cnv_hg19__seg.Level_3'
		file <- sprintf('%s.snp__genome_wide_snp_6__broad_mit_edu__Level_3__segmented_scna_minus_germline_cnv_hg19__seg.seg.txt', study)
	} else if (type == 'meth450') {
		type.str <- 'Merge_methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.Level_3'
		file <- sprintf('%s.methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.data.txt.meth.txt', study)		
	} else if (type == 'meth450genes') {
		type.str <- 'Merge_methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.Level_3'
		file <- sprintf('%s.methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.data.txt.genes.txt', study)		
	} else if (type == 'seg') {
		type.str <- 'Merge_snp__genome_wide_snp_6__broad_mit_edu__Level_3__segmented_scna_minus_germline_cnv_hg19__seg.Level_3'
		file <- sprintf('%s.snp__genome_wide_snp_6__broad_mit_edu__Level_3__segmented_scna_minus_germline_cnv_hg19__seg.seg.txt', study)
	}

	filename <- sprintf("Firehose/stddata__%d_%.2d_%.2d/%s/%d%.2d%.2d/gdac.broadinstitute.org_%s.%s.%d%.2d%.2d00.0.0/%s", 
		year, month, day,study, year, month, day, study, type.str, year, month, day, file)
	return (filename)
}

fhu.filename.analysis <- function(study, type, version = analysis.current.version) {
  day <- as.numeric(version$day)
  month <- as.numeric(version$month)
  year <- as.numeric(version$year)
  if (type == 'cna') {
    type.str <- 'TP.CopyNumber_Gistic2.Level_4'
    file <- 'all_lesions.conf_99.txt'
    file <- 'all_data_by_genes.txt'
  } 
  filename <- sprintf("Firehose/analyses__%d_%.2d_%.2d/%s/%d%.2d%.2d/gdac.broadinstitute.org_%s-%s.%d%.2d%.2d00.0.0/%s", 
                      year, month, day,study, year, month, day, study, type.str, year, month, day, file)
  return (filename)
}


fhu.read.analysis = function(study, type, version = analysis.current.version){
  filename = paths.data(fhu.filename.analysis(study, type, version))
  if (!file.exists(filename)){
    return(NA)
  }
  
  if(type == "cna"){
    cna.data = read.table(filename, sep="\t", header=T, row.names=1, stringsAsFactors = F, check.names = F)
    cna.data = cna.data[,-c(1,2)]
#   cna.data = cna.data[grep("CN values", rownames(cna.data)),]
#   cna.data = cna.data[,-(ncol(cna.data))]
#   colnames(cna.data) = unlist(lapply(colnames(cna.data), function(x){paste(strsplit(x, "\\.")[[1]], collapse="-")}))
#   
#   cna.mat = cna.data[,9:ncol(cna.data)]
#   cna.mat = cbind(unlist(lapply(cna.data$`Wide-Peak-Limits`, function(x){strsplit(x, "\\(")[[1]][1]})), cna.mat)
#   colnames(cna.mat)[1] = "wide_peak_limits"
    
    return(list(analysis_cna = cna.data))
  }
}

fhu.read.analysis.cna = function (study, version = analysis.current.version) {
	return (fhu.read.analysis(study, 'cna', version))
}


fhu.read.mutraw <- function(study, version = fhu.current.version) {
	filename = paths.data(fhu.filename(study, 'mutraw'))
	if (!file.exists(filename))
		return (NA)
	maf.data <- read.delim(filename, header=T, stringsAsFactors=F, sep='\t', check.names = F)
	maf.data = subset(maf.data, Variant_Classification != 'Silent')
	maf.data = subset(maf.data, Tumor_Sample_Barcode != "")
	
 	return (list(mutraw = maf2mat(maf.data), mafraw = maf.data))
}

fhu.read.mut <- function(study, version = fhu.current.version) {
	filename = paths.data(fhu.filename(study, 'mut'))
	if (!file.exists(filename))
		return (NA)
	maf.data <- read.delim(filename, header=T, stringsAsFactors=F, sep='\t', check.names = F)
	maf.data = maf.data[!(maf.data %in% 'Silent'), ]

 	return (list(mut = maf2mat(maf.data)))
}

fhu.read.mutall <- function (study, version = fhu.current.version) {
	mutraw = fhu.read.mutraw(study, version)
	mut = fhu.read.mut(study, version)
	if (is.na(mutraw))
		return (mut)
	mr = tu.tumors(mutraw$mutraw, simplify.barcode = T)
	m  = tu.tumors(mut$mut      , simplify.barcode = T)

	g = intersect(rownames(mr), rownames(m))
	p = setdiff(colnames(mr), colnames(m))
	mutall = cbind(m[g, ], mr[g, p])
 	return (list(mutall = mutall))
}

fhu.read.meth450 <- function(study, version = fhu.current.version) {
	exist1 = file.exists(paths.data(fhu.filename(study, 'meth450')))
	exist2 = file.exists(paths.data(fhu.filename(study, 'meth450genes')))
	if (!exist1 | !exist2)
		return (list (meth450 = NA, meth450genes = NA))

	headers = read.delim(paths.data(fhu.filename(study, 'meth450')), header=T, stringsAsFactors=F, sep='\t', check.names = F, row.names=1, nrows=1)
	classes = c('character', rep('numeric', length(headers)))
	meth450.data <- read.delim(paths.data(fhu.filename(study, 'meth450')), header=T, stringsAsFactors=F, sep='\t', check.names = F, row.names=1, colClasses=classes)
	meth450.genes <- read.delim(paths.data(fhu.filename(study, 'meth450genes')), header=T, stringsAsFactors=F, sep='\t', check.names = F, row.names = 1)
 	colnames(meth450.genes) = c('gene', 'chr')
 	return (list(meth450 = meth450.data, meth450genes = meth450.genes))
}

fhu.read.rppa <- function (study, version = fhu.current.version) {
	filename = paths.data(fhu.filename(study,'rppa'))
	if(!file.exists(filename))
		return (NA)
	rppa.data <- read.delim(filename, row.names=1, header=T, stringsAsFactors=F, sep='\t', check.names = F)
	rownames(rppa.data) <- gsub('.','-',rownames(rppa.data), fixed = T)

	return (list(rppa = rppa.data))
}

fhu.read.clin <- function(study, version = fhu.current.version, followup=NULL) {
	clin.data <- read.delim(paths.data(fhu.filename(study, 'clin')), row.names=1, header=T, stringsAsFactors=F, sep='\t', check.names = F)

#	colnames(clin.data) <- toupper(clin.data['patient.bcrpatientbarcode',])
	colnames(clin.data) <- toupper(clin.data['patient.bcr_patient_barcode',])
	#colnames(clin.data) <- toupper(colnames(data.clin))
	clin.build = clinu.build(clin.data)
	return (list(clin = clin.build))
}

fhu.read.rna <- function (study, version = fhu.current.version) {
	#load GA and Hiseq samples
	print(sprintf('loading %s rnaseqv2-ga...', study))
	ga.num <- NA
	hi.num <- NA
	ct.num <- NA

	ga.filename <- paths.data(fhu.filename(study, 'rna.ga'))
	hi.filename <- paths.data(fhu.filename(study, 'rna.hi'))
	ct.filename <- paths.data(fhu.filename(study, 'rna.counts'))

	if (file.exists(ga.filename)) {
		ga <- try (read.delim(ga.filename,header=T, check.names=F, row.names=1, sep='\t', stringsAsFactors=F, fill=T))
		ga <- ga[-1, ]
		ga.num <- matrix(as.numeric(as.matrix(ga)), nrow=nrow(ga), dimnames = list(rownames(ga), colnames(ga)))
	} else {
		print(sprintf('could not find GA rna files for %s', study))
	}

	print(sprintf('loading %s rnaseqv2-hiseq...', study))
	if (file.exists(hi.filename)) {
		hi <- read.delim(hi.filename,header=T, check.names=F, row.names=1, sep='\t', stringsAsFactors=F, fill=T)
		hi <- hi[-1, ]
		hi.num <- matrix(as.numeric(as.matrix(hi)), nrow=nrow(hi), dimnames = list(rownames(hi), colnames(hi)))
	} else {
		print(sprintf('could not find HISeq rna files for %s', study))
	}

	print(sprintf('loading %s rnaseqv2 counts...', study))
	if (file.exists(ct.filename)) {
		ct <- read.delim(ct.filename,header=T, check.names=F, row.names=1, sep='\t', stringsAsFactors=F, fill=T)
		ct <- ct[-1, ]
		ct.num <- matrix(as.numeric(as.matrix(ct)), nrow=nrow(ct), dimnames = list(rownames(ct), colnames(ct)))

		genes <- tu.genes(ct.num)
		ct.num <- ct.num[!duplicated(genes), ]
		rownames(ct.num) <- tu.genes(ct.num)
		ct.num = round(ct.num)
	} else {
		print(sprintf('could not find CountV2 rna files for %s', study))
	}

	return (list(ga = ga.num, hi = hi.num, ct = ct.num))
}

fhu.read.cna <- function (study,  version = fhu.current.version) {
	cna.data <- read.delim(paths.data(fhu.filename(study, 'cna')),  header=T, stringsAsFactors=F, sep='\t', check.names = F)
	x = cnu.CNmatrix(cna.data, load.ref=T)
	return (list(cna = x))
}

fhu.read.cnag <- function (study,  version = fhu.current.version) {
	cna.data <- read.delim(paths.data(fhu.filename(study, 'cnag')),  header=T, stringsAsFactors=F, sep='\t', check.names = F)
	x = cnu.CNmatrix(cna.data, load.ref=T)
	return (list(cnag = x))
}


fhu.read.seg <- function (study, version = fhu.current.version) {
	filename = paths.data(fhu.filename(study, 'seg', version = fhu.current.version))
	seg = read.table(filename, stringsAsFactors = F, header = T, sep = '\t', check.names = F) 

	## turn segmentation file into a by-gene matrix
	seg.tum = seg[tu.istumor(seg$Sample), ]
	segmat_by_gene = cnu.CNmatrix(seg.tum, load.ref = T)

	return (list (seg = seg, segmat_by_gene = segmat_by_gene))
}


fhu.update.study <- function (study, version = fhu.current.version, skip.load = FALSE, update = NA) {
	# load all fh data, and save a structure for quick upload
	ptm <- proc.time()

	day <- as.numeric(version$day)
	month <- as.numeric(version$month)
	year <- as.numeric(version$year)
	filename <- paths.processed(sprintf('firehose_%.4d_%.2d_%.2d_%s.rds', year, month, day, study))
	print(filename)

	if (!skip.load) {
		if (file.exists(filename)) {
			print(sprintf('%s, fhu.get.study: loading existing file %s', study, filename))
			tcga = readRDS(filename)
		}
	}

	if (!is.na(update)) {
		datas = update
	} else {
		datas = list (	list(name = 'RNAseqV2',	fun = fhu.read.rna),
						list(name = 'Mutation',	fun = fhu.read.mut),
						list(name = 'MutationRaw',	fun = fhu.read.mutraw),
						list(name = 'MutationAll',	fun = fhu.read.mutall),
						list(name = 'RPPA',		fun = fhu.read.rppa),
						list(name = 'Clinical',	fun = fhu.read.clin),
						list(name = 'Meth',		fun = fhu.read.meth450),
						list(name = 'Seg', 		fun = fhu.read.seg),
						list(name = 'Copy num germ', fun = fhu.read.cnag),
						list(name = 'Analysis_cna', fun = fhu.read.analysis.cna ))
		tcga <- list()

	}


	for (data in datas) {
		print(sprintf('%s, loading %s data...', study, data$name))
		dat  <- tryCatch(data$fun(study), error=errorFileNotFoundWarning)
		if (!is.list(dat)) {
			print(sprintf('%s: error: %s not found', study, data$name))
		} else {
			for (field in names(dat)) {
				tcga[[field]] = dat[[field]]
			}
		}
	}

	print(sprintf('%s: fhu.get.study: saving file %s', study, filename))
	saveRDS(tcga, file = filename)

	print(proc.time() - ptm)
	return (NA)
}



fhu.read.study <- function (study, version = fhu.current.version) {
	# load all fh data, and save a structure for quick upload
	ptm <- proc.time()

	day <- as.numeric(version$day)
	month <- as.numeric(version$month)
	year <- as.numeric(version$year)
	filename <- paths.processed(sprintf('firehose_%.4d_%.2d_%.2d_%s.rds', year, month, day, study))
	print(filename)

	if (file.exists(filename)) {
		print(sprintf('fhu.get.study: loading existing file %s', filename))
		return(readRDS(filename))
	} else {
		return (NA)
	}
}

fhu.report <- function (studies = fhu.studies(), mc.cores=1) {
	datasets = c('ct', 'ga', 'hi', 'mut', 'cnag', 'rppa', 'clin', 'meth450', 'meth450genes', 'fields')
	#report = matrix("", ncol=length(datasets), nrow = length(studies))
	#colnames(report) = datasets
	#rownames(report) = studies
	reports = mclapply (studies, 
		function (s) {
			report = rep("", length(datasets))
			names(report) = datasets
			print(sprintf('processing %s', s))
			study = fhu.read.study(s)
			fields = names(study)
			report['fields'] = paste(fields, sep='', collapse=', ')
			reportfields = intersect(fields, datasets)
			for (f in reportfields) {
				if (all(is.na (study[[f]]))) {
					report[f] = "NA"
				}
				else {
					if (f == 'clin')
						r = sprintf('%d', length(study[[f]]$patient$bcr_patient_barcode))
					else
						r = sprintf('(%d,%d)', nrow(study[[f]]), ncol(study[[f]]))
					report[f] = r
				}
			}; 
			return(report); 
		}, 
		mc.cores = mc.cores)

	reports = do.call(rbind, reports)
	rownames(reports) = studies
	return (reports)
}


fhu.plot <- function (study.df, x, y, muts = NA, xlog = F, ylog = F, grep.xy.fixed = T,  grep.muts.fixed = T) {
# example: fhu.plot(tcga[['BRCA']], 'COL1A1','FN1_Fibronectin_p', c('PIK3CA_m', 'PTEN_m','TP53_m'), grep.xy.fixed = F, grep.muts.fixed = T, xlog=T, ylog=F) 

	p = Reduce(intersect, list(colnames(study.df$rsem), colnames(study.df$mut), colnames(study.df$rppa)))
	m = study.df$mut[,p]
	rownames(m) <- paste(rownames(m),'m', sep = '_')
	r = rbind(study.df$rsem[,p], m, study.df$rppa[,p])
#	r = rbind(study.df$rsem[,p], study.df$mut[,p], study.df$rppa[,p])
#	p = Reduce(intersect, list(colnames(study.df$rsem), colnames(study.df$mut)))
#	r = rbind(study.df$rsem[,p], study.df$mut[,p])

	x.ix = grep (x, rownames(r), fixed = grep.xy.fixed)[1]
	y.ix = grep (y, rownames(r), fixed = grep.xy.fixed)[1]

	fx <- function (x) {return (x)}
	fy <- function (x) {return (x)}
	if (xlog) fx <- function (x) {log2(1+x)}
	if (ylog) fy <- function (x) {log2(1+x)}

	leg.names = c()
	colors = c('red', 'blue', 'green', 'yellow', 'purple', 'black')
	plot (fx(as.numeric(r[x.ix, ])), fy(as.numeric(r[y.ix,])), pch=20, cex= 0.3, xlab = x, ylab = y)
	for (i in 1:length(muts)) {
		mix = grep (muts[i], rownames(r), fixed = grep.muts.fixed)
		mix.pat = colnames(r)[r[mix,] == 1]
		points (fx(as.numeric(r[x.ix,mix.pat])), fy(as.numeric(r[y.ix,mix.pat])), col = colors[i], cex = 0.5+i/3)
		leg.names = c(leg.names, muts[i])
	}
	legend('topleft', leg.names, col = colors[1:length(muts)], cex = 1, pch= 1)
}
