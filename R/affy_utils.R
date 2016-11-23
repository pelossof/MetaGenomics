# Set of functions for affy data loading/normalization
# Manu N Setty
# 11/27/2009
# revides by Rafi Pelossof 5/2015

# Revised         Comments
# 04/01/2010      Added probes.to.symbols function
# 02/10/2011      Speed up probes.to.symbols by calling aafSymbol only once
# Function to include all the required libraries
#

include.affy.libraries <- function () {
  require (affy)
  require (annaffy)
  require (limma)
}

# Function to read the CEL files
# fnames   Vector of file names to be loaded
# Returns  the object returned by ReadAffy function
read.cel.files <- function (fnames) {

  include.affy.libraries ()
  
  time.start <- get.time ()
  cel <- ReadAffy(filenames=fnames)
  time.end <- get.time()
  show(sprintf("Minutes to read in cel files: %0.1f", (time.end - time.start)/60))
  return(cel)
}


# Function to normalize affy batch data.
# Uses expresso
# RMA background correction,
# Quantile normalization
# Median polish summarisation
# Uses only Perfect Match probes

# cel          ReadAffy batch object
# colnames     (Optional). Column name specifier for resulting matrix.
# rownames     (Optional). Row name specifier for resulting matrix.
#              Can also be the name of the chip. If chip  is specified, the
#              corresponding row gene annotations are used for median
#              summarization

# Returns normalized data matrix: No. of genes x No. of samples
normalize.cel.data <- function (cel, colnames=NA, rownames=NA) {
  include.affy.libraries ()
  
  time.start <- get.time ()
  normExprs <- expresso(cel, bgcorrect.method = 'rma',
                        normalize.method = 'quantiles',
                        pmcorrect.method = 'pmonly',
                        summary.method = 'medianpolish')

  time.end <- get.time ()
  show(sprintf("Minutes to normalize cel files: %0.1f", (time.end - time.start)/60))
  
  # Make a temporary object
  save (normExprs, file="__exprs_data.Rdata")
  exprs.data <- exprs (normExprs)
  
  if (hasArg (colnames))
    colnames (exprs.data) <- colnames

  if (hasArg (rownames)) {
    # Check if the rownames is a chip name
    if (length (rownames)==1) {
      time.start <- get.time ()

      gene.names <- probes.to.symbols (rownames (exprs.data), rownames)
      aggregate.data <- aggregate (exprs.data, by=list (gene.names), FUN=median)
      #duplicateCorrelation solves the multiprobe problem
      #fit <- lmFit(MA, design, ndups = 2, spacing = 1, cor = corfit$consensus)
      rownames (aggregate.data) <- aggregate.data[[1]]
      aggregate.data <- aggregate.data[,-1]
      exprs.data <- round (aggregate.data, digits = 4)

      time.end <- get.time ()
      show(sprintf("Minutes to summarized data from annotations: %0.1f", (time.end - time.start)/60))
    }
    else
      rownames (exprs.data) <- rownames
  }

  return (exprs.data)
}
  

probes.to.symbols <- function (probes, chip) {
  include.affy.libraries ()
  symbols <- sapply (aafSymbol (probes, chip), function (x) { ifelse (length (x!=0), x, NA) } )
  names (symbols) <- probes
  return (symbols)
}


get.time <- function () 
  as.numeric(format(Sys.time(), "%s"))

##########################################
### Added by Rafi Pelossof Sep 26 2013 ###
##########################################

hist.gene <- function (x, y, test.gene, gene.rownames = NULL) {
  if(!hasArg("gene.rownames"))
    gene.rownames <- rownames(x)
  x1 <- as.vector(x[grepl(test.gene,gene.rownames), y==1])
  x0 <- as.vector(x[grepl(test.gene,gene.rownames), y==0])
  breaks <- seq(floor(min(c(x1,x0))),ceiling(max(c(x1,x0))), length=20)
  p.mut <- hist(x[grepl(test.gene,gene.rownames), y==1],breaks = breaks)
  p.wt  <- hist(x[grepl(test.gene,gene.rownames), y==0],breaks = breaks)
  p.mut$counts <- p.mut$counts/sum(p.mut$counts)
  p.wt$counts <- p.wt$counts/sum(p.wt$counts)
  x.min <- min(c(p.mut$breaks, p.wt$breaks))
  x.max <- max(c(p.mut$breaks, p.wt$breaks))
  y.max <- max(c(p.mut$counts, p.wt$counts))*1.1
  plot( p.mut, col=rgb(1,0,0,1/4), xlim=c(x.min,x.max), ylim=c(0,y.max))
title(sub=test.gene)
  plot( p.wt, col=rgb(0,0,1,1/4), xlim=c(x.min,x.max), ylim=c(y.max), add=TRUE)

  return (ks.test(x0, x1))
}

ks <- function(x, y, test.gene, gene.rownames) {
  x1 <- as.vector(x[grepl(test.gene,gene.rownames), y==1])
  x0 <- as.vector(x[grepl(test.gene,gene.rownames), y==0])
  #print(paste(length(x0), "/", length(x1)))
  return (ks.test(x0, x1))

}
