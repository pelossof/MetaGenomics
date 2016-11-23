library(org.Hs.eg.db)
library(GenomicFeatures)

source('~/paths.R')
source.util('tcga')
source.util('maf')

# Function to replicate the CN.summary track in IGV when a seg file is loaded
# Returns a GRanges object
#    segFile = matrix with seven columns: (sample, chrom, loc.start, loc.end, num.mark, num.informative, seg.mean)
cnu.IGVsummary <- function(segFile, skip.load=F){
  library(GenomicRanges)
  library(biomaRt)
  
  # Get chromosome lengths from UCSC
  if(skip.load){
    chromLengths <- load(paths.data('processed/cnu.chrLen.rds'))
  }else{
    con <- gzcon(url(paste("http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/chromInfo.txt.gz", sep="")))
    txt <- readLines(con)
    chroms <- read.table(textConnection(txt))
    chroms <- chroms[grep("_", chroms[,1], invert=T), 1:2]
    chroms <- chroms[order(as.numeric(gsub("chr","",chroms[,1]))),]
    chroms <- chroms[chroms[,1] != "chrM",]
    chromLengths <- chroms[,2]
  }
  
  binSize <- 200000
  threshold <- 0.1
  results <- c()
  
  numSamples <- length(unique(segFile[,'Sample']))
  
  # for each chromosome
  for (i in 1:length(chromLengths)){
    chr <- i
    chrLen <- chromLengths[i]
    bins <- seq(1,chrLen, binSize)
    bins <- cbind(bins,bins+binSize-1, 0,0)
    colnames(bins) <- c("binStart", "binEnd", "amp","del")
    bins[nrow(bins),2] <- chrLen
    
    segs <- segFile[segFile$chr == chr,]
    
    # For each segment on that chromosome: 
    if (nrow(segs) != 0){
      for (j in 1:nrow(segs)){
        segStart <- segs$loc.start[j]
        segEnd <- segs$loc.end[j]
        
        startBin <- ceiling(segStart / binSize)
        endBin <- ceiling(segEnd / binSize)
        
        if(endBin > nrow(bins)){
          endBin <- nrow(bins)
        }
        
        # For each bin spanned by that segment:
        if (startBin <= endBin){
          for (bin in seq(startBin, endBin)){
            binStart <- bins[bin,1]
            binEnd <- bins[bin,2]
            
            weight <- 1
            if (segEnd < binEnd){
              weight <- (segEnd - max(segStart, binStart)) / binSize
            }else if (segStart > binStart){
              weight <- (min(segEnd, binEnd) - segStart) / binSize
            }
            
            if(is.finite(segs$seg.mean[j])){
              if (segs$seg.mean[j] > threshold){
                bins[bin, 3] <- bins[bin, 3] + weight
              }else if (segs$seg.mean[j] < -1*threshold){
                bins[bin, 4] <- bins[bin, 4] - weight
              }
            }
          }
        }
      }
    }
    results <- rbind(results, cbind(chroms[chr,1], bins))
  }
  
  segments.GR <- GRanges(seqnames <- results[,1], ranges <- IRanges(as.numeric(results[,2]), end=as.numeric(results[,3])), strand='*',
                        amp <- as.numeric(results[,4])/numSamples, del=as.numeric(results[,5])/numSamples)
  return(segments.GR)
}


# Makes gene reference object needed for cnu.CNmatrix
cnu.makeRef <- function(load.ref = F){
  filename = paths.processed('cnu.tx_hg19.rds')
  if (file.exists(filename) & load.ref) {
      print('cnu.makeRef: loading CN reference...')
      ref <- readRDS(filename)
      return(ref)
  } else {
    print('cnu.makeRef: creating CN reference...')
    hg19.ts <- makeTranscriptDbFromUCSC(genome <- 'hg19', tablename <- 'refGene')
    cds = unlist(exonsBy(hg19.ts, 'tx', use.names=T))
    cds$tx_id <- names(cds)
    geneName <- select(org.Hs.eg.db, keys=names(cds), columns=c("SYMBOL"), keytype="REFSEQ")
    geneName[is.na(geneName$SYMBOL),2] <- cds$tx_id[is.na(geneName$SYMBOL)]
    cds$symbol <- geneName$SYMBOL
    ref <- unique(cbind(cds$symbol, as.character(seqnames(cds)), as.character(strand(cds))))
    ref <- ref[grep("_", ref[,2], invert=T),]
    rownames(ref) <- ref[,1]
    
    start_pos <- aggregate(start(cds), by=list(cds$symbol), min)
    end_pos <- aggregate(end(cds), by=list(cds$symbol), max)
    rownames(start_pos) <- start_pos[,1]
    rownames(end_pos) <- end_pos[,1]
    
    ref = cbind(ref, NA, NA)
    colnames(ref) <- c("ID", "chr","strand","start","end")
    ref[, "start"] <- start_pos[rownames(ref),2]
    ref[, "end"] <- end_pos[rownames(ref),2]
    
    gene_chr = unique(cbind(cds$symbol, as.character(seqnames(cds)), as.character(strand(cds))))
    bad_genes = names(which(table(gene_chr[,1]) > 1))
    bad_genes = c("NR_120535", bad_genes)
    ref = ref[!ref[,"ID"] %in% bad_genes,]

    ref.gr = GRanges(seqnames=Rle(ref[,"chr"]),
                     ranges=IRanges(as.numeric(ref[,"start"]), end=as.numeric(ref[,"end"]), names=ref[,"ID"]),
                     strand=ref[,"strand"])
    
    
    saveRDS(ref.gr, file = filename)
    return(ref.gr)
  }
}

# Function to create a table of seg_means for patient vs gene/miR
#    segFile = matrix with seven columns: (sample, chrom, loc.start, loc.end, num.mark, num.informative, seg.mean)
#    geneLookup = gene/miR/probe reference: matrix with 5 columns: (id, chromosome, strand, start, end), from cnu.makeRef
#
# expecting segFile data frame with the columns: Sample, Chromosome, Start, End, Segment_Mean, each row is a CNA event
# geneRef should be a ...
cnu.CNmatrix <- function(segFile, geneRef = NA, progress=T, load.ref=F){
  if (is.na(geneRef)) geneRef <- cnu.makeRef(load.ref=load.ref)
  
  if (length(grep("chr", segFile[,'Chromosome'])) == 0){
    segFile[,'Chromosome'] <- paste("chr", segFile[,'Chromosome'], sep="")
  }
  
  segFile[,"Chromosome"] = gsub("chr23", "chrX", segFile[,"Chromosome"])
  segFile[,"Chromosome"] = gsub("chr24", "chrY", segFile[,"Chromosome"])
  
  sampleNames <- unique(segFile[,'Sample'])
  gene_patient <- matrix(0, nrow=length(sampleNames), ncol=length(geneRef))
  rownames(gene_patient) = sampleNames
  colnames(gene_patient) = names(geneRef)
  if(progress){
    pb <- txtProgressBar(min = 1, max = length(sampleNames), width=100, style=3)
  }
  for (i in 1:length(sampleNames)){
    sample <- sampleNames[i]
    #sampleMat <- as.matrix(segFile[segFile[,'Sample'] == sample,])  # subset of segFile for this sample
    sampleMat <- segFile[segFile[,'Sample'] == sample,]  # subset of segFile for this sample
    
    sample.gr = GRanges(seqnames=Rle(sampleMat[,"Chromosome"]),
                        ranges=IRanges(as.numeric(sampleMat[,"Start"]), end=as.numeric(sampleMat[,"End"])),
                        segment_mean=sampleMat[,"Segment_Mean"])
    hits = as.data.frame(findOverlaps(geneRef, sample.gr))
    gene_patient[i,hits[!hits[,1] %in% hits[duplicated(hits[,1]),1],1]] = as.numeric(sample.gr$segment_mean[hits[!hits[,1] %in% hits[duplicated(hits[,1]),1],2]])
    
    if(progress){
      setTxtProgressBar(pb, i)
    }
  }
  if(progress) print("")
  
  return(tu.tumors(t(gene_patient)))
}

cnu.toGene = function(regions, geneRef = NA, load.ref=F){
  if (is.na(geneRef)) geneRef <- cnu.makeRef(load.ref=load.ref)
  
  genes = rep(list(list()), length(regions))
  names(genes) = regions

  for(region in regions){
    chr = strsplit(region, ":")[[1]][1]
    start = as.numeric(strsplit(region, "[-:]")[[1]][2])
    end = as.numeric(strsplit(region, "-")[[1]][2])
    region.genes = names(geneRef)[as.character(seqnames(geneRef)) == chr & start(geneRef) < end & end(geneRef) > start]
    
    genes[[region]] = region.genes
  }

  return(genes)  
}


diff.cnag = function(cohorts, study, load.ref=F){
  dat = study$cna
  regions = dat[,1]
  dat = tu.tumors(dat[,-1], simplify.barcode = 1)
#  colnames(dat) = unlist(lapply(colnames(dat), function(x){paste(strsplit(x, "-")[[1]][1:3], collapse="-")}))
  
  matrix1 = dat[,cohorts[[1]][cohorts[[1]] %in% colnames(dat)]]
  matrix2 = dat[,cohorts[[2]][cohorts[[2]] %in% colnames(dat)]]
  
  results = list()
  for(i in 1:nrow(dat)){
    results[[i]] = t.test(as.numeric(matrix1[i,]), as.numeric(matrix2[i,]))
  }
  pval = unlist(lapply(results, function(x){x$p.value}))
  padj = p.adjust(pval)
  mean1 = unlist(lapply(results, function(x){x$estimate[1]}))
  mean2 = unlist(lapply(results, function(x){x$estimate[2]}))
  
  res.matrix = cbind(as.character(regions), substr(rownames(dat), 1, 3), mean1, mean2, pval, padj)
  colnames(res.matrix) = c("Region","CN", "MeanX", "MeanY", "pvalue", "padj")
  rownames(res.matrix) = NULL
  res.matrix = res.matrix[order(res.matrix[,'padj']), ]
  return(as.data.frame(res.matrix))
}






