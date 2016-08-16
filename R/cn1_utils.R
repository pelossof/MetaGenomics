

# http://www.psytky03.com/2015/04/analysis-of-agilent-array-cgh-data.html

# The code is for the analysis of Agilent high density Array-CGH data 
# to identify copy number variation regions. 
# A Bioconductor package DNAcopy is used. 
# The array analyzed is SurePrint G3 Mouse CGH Microarray 1x1M_1M, 
# but the code should work as well for other arrays

require("DNAcopy")

#set the working directory where the raw data saved


# define a function "convert_rawdata"
# This function is used to extract raw logR ratio data and chromosome positions to a data frame subject
cnu.loadArray <- function (raw_file)
{

	#raw_file = "/Volumes/Pegasus24/data/JGA/CGH/R01_Biopsy_Before_Treatment/ID/107-B.txt"
	print(sprintf('loading array %s', raw_file))
	load.error = F

	temp.file = tryCatch(read.table(file = raw_file, skip = 9,  header=T, stringsAsFactors =FALSE, fill = T, sep = '\t'), 
		error = function(cond){print(sprintf('Error loading %s, condition: %s', raw_file, cond)); load.error = T; return (NA)})

	if (class(temp.file) != 'data.frame')
		return (NA)

	print(sprintf('loaded array %s, nrows = %d', raw_file, nrow(temp.file)))


	valid.probes = grepl ('chr[0-9]*[XY]*:', temp.file$SystematicName)
	valid.arr = temp.file[valid.probes, ]

	temp.pos = as.vector(sapply(strsplit(valid.arr$SystematicName,":"), "[", 2))
	final.pos = as.numeric(sapply(strsplit(temp.pos,"-"),"[",1))
	temp.chr = as.vector(sapply(strsplit(as.character(valid.arr$SystematicName),":"), "[", 1))
	final.chr <- gsub("chr", "", temp.chr)
	final.chr[final.chr == 'X'] = 23
	final.chr[final.chr == 'Y'] = 24
	final.chr = as.numeric(final.chr)

	# Prepare the data frame 

	obj_df <- data.frame(logR=valid.arr$LogRatio * log10(2),
	                   Chr=final.chr,
	                   Pos=as.numeric(final.pos), stringsAsFactors = F)
	obj_df = obj_df[with(obj_df, order (Chr, Pos)), ]
	return(obj_df)
}

cnu.arr2obj = function (arr) {
	CNA(arr$logR, arr$Chr, as.numeric(arr$Pos), data.type="logratio", sampleid=1)
}

cnu.smoothObj = function (obj, ...) {
	return (segment(smooth.CNA(obj),...))
}

cnu.smoothObj2df = function(obj) {
	segment.smoothed.CNA.object = obj
	
	# Save the output data 
	segments <- segment.smoothed.CNA.object$output
#	signal_data <- as.data.frame(segment.smoothed.CNA.object$data)

	# reset the positions
	#colnames(segments)
	segments$loc.start = segments$loc.start*1000
	segments$loc.end = segments$loc.end*1000
#	signal_data$maploc = signal_data$maploc*1000

	# Log10 -> Log2 convertion 
	# Log2 is generally used to interpret the data, 0.5 indicates 3 copy duplication.
#	signal_data[,c(3,4)] <- log2(10)*signal_data[,c(3,4)]
	segments$seg.mean <-  log2(10)*segments$seg.mean

	return (segments)	

	# Save the segmentation file
	#write.csv(segments, file="segmentation.csv")
}

cnu.fname2seg = function (filename) {
	arr = cnu.loadArray(filename)
	return (cnu.arr2seg(arr))
}

cnu.arr2seg = function (arr, name = 'sample') {
	obj = cnu.arr2obj(arr)
	obj.smooth = cnu.smoothObj(obj, alpha = 0.001)
	segments = cnu.smoothObj2df(obj.smooth)
	segments$ID = name
	return (segments)
}

cnu.plots = function (obj.smooth) {
	plot(segment.smoothed.CNA.object,plot.type="s")
	plot(segment.smoothed.CNA.object, plot.type="p")
	plot(segment.smoothed.CNA.object, plot.type="c")
	plot(segment.smoothed.CNA.object, plot.type="w")
}

