options(stringsAsFactors = F)
library(GenomicFeatures)
library(geneModels)
library(org.Hs.eg.db)
library(Rsamtools)
library(multicore)
library(annotate)
library(data.table)

####################################################################################################
#This function extracts the annotation for the gene and exons to be used for getting deseq gene counts 
#and dexseq read counts.
#It takes the flattened genome as the input
#It returns a list with two objects: 
#1. gene.annotation  <- GRangesList of exons by gene 
#2. exon.annotation <- GRanges of exons
####################################################################################################
ExtractAnnotation <- function (annotation.object.location) {
  
  gr.annotation <- readRDS(annotation.object.location)
  
  #Remove the ranges coressponding to intron, intergeneic ,overlap, extended 3'UTR and extened 5'UTR
  gr.annotation <- gr.annotation[!(gr.annotation$exon.anno %in% c("intron", "intergenic", "overlap", "utr5*", "utr3*"))]
  
  #Reducing the granges based on entrez gene id(using annotated ranges) and converting it back to granges object
  ar.annotation <- geneModels:::AnnotatedGRanges(gr.annotation, as.character(gr.annotation$entrez.id))
  ar.annotation <- reduce(ar.annotation)
  gr.annotation <- as(ar.annotation, "GRanges")
  gr.annotation <- gr.annotation[order(gr.annotation$Symbol)]
  
  #Maps the entrez ids to gene symbols
  vec.symbol <- unlist(lookUp(values(gr.annotation)$Symbol, "org.Hs.eg", "SYMBOL"), use.names = F)
  #metadata frame with entrez id and symbol columns
  df.metadata <- data.frame(entrez.id = values(gr.annotation)$Symbol, symbol = vec.symbol, stringsAsFactors = FALSE)
    
  #getting the frequency of each entrez id(unfortunately it gets sorted)
  df.num.genes <- as.data.frame(table(df.metadata$entrez.id), responseName = "freq")
  colnames(df.num.genes) <- c("entrez.id", "freq")
  
  #Keeping the same order of entrez id as in original
  dt.num.genes <- data.table(df.num.genes, key = c("entrez.id"))
  df.num.genes <- data.frame(dt.num.genes[unique(df.metadata$entrez.id), ])
  
  #getting a new column that defines the exon
  num.genes <- nrow(df.num.genes)
  vec.exon.names <- paste("exon", unlist(mapply(seq, from = rep(1, times = num.genes) , to = df.num.genes$freq), use.names = F), 
                          sep = "")
  df.metadata <- cbind(df.metadata, exon = vec.exon.names)
  
  #Assigning the newly defined metadata to the granges object
  values(gr.annotation) <- df.metadata
  
  #Converting the granges object to granges list by using the entrez id as the split factor
  vec.exon.times <- unlist(mapply(rep, x = 1:num.genes, each = df.num.genes$freq), use.names = F)
  grl.annotation <- split(x = gr.annotation, f = vec.exon.times)
  names(grl.annotation) <- unique(df.metadata$entrez.id)
  
  return(list(gene.annotation = grl.annotation, exon.annotation = gr.annotation))
}


####################################################################################################
## This function takes the input as the path of the bam file for single end reads, gene annotation,  
## exon annotation and number of cores. It returns a list with the reads counts for each gene and exon
## ARGUMENTS
##1. file.path <- path for bam file
##2. grl.gene.annotation <- GRangesList that defines the gene annotations
##3. grl.exon.annotation <- GRanges that defines the exon annotations
##4. num.cores <- the number of cores for parallel processing
##
## VALUES
##1. list of gene.counts and exon.counts
####################################################################################################
ReadBamSingleEndUnstranded <- function(file.path, grl.gene.annotation, gr.exon.annotation, num.cores) {
  
  cat("Reading File", file.path, "\n")
  
  vec.chr <- c("chrX", "chrY", "chrM", paste("chr", seq(1, 22), sep = ""))
  ls.counts <- mclapply(seq(1,length(vec.chr)), function(x) {
    
    #Specifying range of a particular chromosme
    gr.chr <- GRanges(vec.chr[x], IRanges(start = 1, end = 536870912))
    #Defining the parameters to be read from the bam file
    parameter <- ScanBamParam(tag="NH", which = gr.chr)
    ga.alignment <- readGappedAlignments(file.path, param = parameter, use.names = TRUE)
    
    #finding the width of the reads and dividing by 2
    overlap <- floor(qwidth(ga.alignment[1])/2)
    #finding unique hits
    vec.unique.hit <- values(ga.alignment)$NH == 1L
    
    #Get all the reads with unique alignment 
    ga.unique.alignment <- ga.alignment[vec.unique.hit]
    
    #get the counts for overlapping reads
    vec.gene.counts <- unname(assays(summarizeOverlaps(features = grl.gene.annotation, reads = ga.unique.alignment, 
                                                minoverlap = overlap, ignore.strand = T,
                                                mode="IntersectionNotEmpty"))$counts)
    
    vec.exon.counts <- unname(assays(summarizeOverlaps(features = gr.exon.annotation, reads = ga.unique.alignment, 
                                                minoverlap = overlap, ignore.strand = T,
                                                mode="IntersectionNotEmpty"))$counts)
    
    return(list(gene = vec.gene.counts, exon = vec.exon.counts))
  }, mc.cores = num.cores, mc.cleanup = T)
  
  #Sums up the counts across all the chromosomes
  vec.gene.counts <- rowSums(do.call(cbind, lapply(seq(1, length(ls.counts)), function(x) ls.counts[[x]]$gene)))
  vec.exon.counts <- rowSums(do.call(cbind, lapply(seq(1, length(ls.counts)), function(x) ls.counts[[x]]$exon)))
  
  return(list(gene.counts = vec.gene.counts, exon.counts = vec.exon.counts))
}


####################################################################################################
## This function takes the input as the path of the bam file for pairedend reads, gene annotation,  
## exon annotation and number of cores. It returns a list with the reads counts for each gene and exon
## ARGUMENTS
##1. file.path <- path for bam file
##2. grl.gene.annotation <- GRangesList that defines the gene annotations
##3. grl.exon.annotation <- GRanges that defines the exon annotations
##4. num.cores <- the number of cores for parallel processing
##
## VALUES
##1. list of gene.counts and exon.counts
####################################################################################################
ReadBamPairedEndUnstranded <- function(file.path, grl.gene.annotation, gr.exon.annotation, num.cores) {
  
  cat("Reading File", file.path, "\n")
  
  vec.chr <- c("chrX", "chrY", "chrM", paste("chr", seq(1, 22), sep = ""))
#  ls.counts <- mclapply(seq(1,length(vec.chr)), function(x) {
 ls.counts <- mclapply(seq(1,length(vec.chr)), function(x) {    
    #Specifying range of a particular chromosme
    gr.chr <- GRanges(vec.chr[x], IRanges(start = 1, end = 536870912))
    #Defining the parameters to be read from the bam file
    parameter <- ScanBamParam(tag="NH", flag= scanBamFlag(isProperPair = T), which = gr.chr, what=c("flag", "mrnm", "mpos"))
    ga.alignment <- readGappedAlignments(file.path, param = parameter, use.names = TRUE)
    
    #finding the width of the reads and dividing by 2
    overlap <- floor(qwidth(ga.alignment[1])/2)
    #finding unique hits
    vec.unique.hit <- values(ga.alignment)$NH == 1L
    
    #Get all the reads with unique alignment and then pair it
    ga.unique.alignment <- ga.alignment[vec.unique.hit]
    gap.alignment <- makeGappedAlignmentPairs(ga.unique.alignment)
    
    #get the counts for overlapping reads
    vec.gene.counts <- unname(assays(summarizeOverlaps(features = grl.gene.annotation, reads = gap.alignment, 
                                                minoverlap = overlap, ignore.strand = T,
                                                mode="IntersectionNotEmpty"))$counts)
    
    vec.exon.counts <- unname(assays(summarizeOverlaps(features = gr.exon.annotation, reads = gap.alignment, 
                                                minoverlap = overlap, ignore.strand = T,
                                                mode="IntersectionNotEmpty"))$counts)
    
    return(list(gene = vec.gene.counts, exon = vec.exon.counts))
  }, mc.cores = num.cores, mc.cleanup = T)
  
  #Sums up the counts across all the chromosomes
  vec.gene.counts <- rowSums(do.call(cbind, lapply(seq(1, length(ls.counts)), function(x) ls.counts[[x]]$gene)))
  vec.exon.counts <- rowSums(do.call(cbind, lapply(seq(1, length(ls.counts)), function(x) ls.counts[[x]]$exon)))
  
  return(list(gene.counts = vec.gene.counts, exon.counts = vec.exon.counts))
}

##########################################################################################################
#This function extracts the name of the sample
#I change the function according to thhe naming convention for the sample
##########################################################################################################
ExtractSampleName <- function(file.name){
  sample.name <- unlist(strsplit(bam.file, "-"), use.names = F)[2]
  final.sample.name <- gsub("_accepted_hits.bam", "", sample.name)
  return(final.sample.name)
}


main <- function () {
#########################################################################################################
#Running accordingly
########################################################################################################
annotation.object.location <- "/ifs/e63data/leslielab/singh/annotation_data/annotated.genome.hg19.rds"
ls.annotations <- ExtractAnnotation(annotation.object.location)
grl.gene.annotation <- ls.annotations$gene.annotation
gr.exon.annotation <- ls.annotations$exon.annotation


m.gene.counts <- NULL
m.exon.counts <- NULL
vec.sample.names <- NULL
num.cores <- 20
file.names = list.files(path = "/ifs/e63data/leslielab/singh/yilong_project/alignments", pattern = "accepted_hits.bam$", full.names = T, recursive = T)

for(bam.file in file.names){
  
  sample.name <- ExtractSampleName(bam.file)
  ls.counts <- ReadBamPairedEndUnstranded (bam.file, grl.gene.annotation, gr.exon.annotation, num.cores)
  
  m.gene.counts <- cbind(m.gene.counts, ls.counts$gene)
  m.exon.counts <- cbind(m.exon.counts, ls.counts$exon)
  vec.sample.names <- c(vec.sample.names, sample.name)
  
}  

colnames(m.gene.counts) <- vec.sample.names
colnames(m.exon.counts) <- vec.sample.names

se.gene <- SummarizedExperiment(rowData = grl.gene.annotation, assays = SimpleList(counts = m.gene.counts))
se.exon <- SummarizedExperiment(rowData = gr.exon.annotation, assays = SimpleList(counts = m.exon.counts))


gene.file.location <- "/ifs/e63data/leslielab/singh/yilong_project/data/counts/se_gene.rds"
exon.file.location <- "/ifs/e63data/leslielab/singh/yilong_project/data/counts/se_exon.rds"
saveRDS(se.gene, file = gene.file.location)
saveRDS(se.exon, file = exon.file.location)

}