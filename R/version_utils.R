#
# Rafi Pelossof, 2015, MSKCC
# version_utils.R
#
#	dump everything to results/versionX/...
#


ver.new = function (name = NA, v = NA) {
	if(is.na(v))
		return(paste(format(Sys.time(), "%Y%m%d_%H%M%S"), ifelse(is.na(name),'',sprintf('_%s', name)), sep = ''))
	else
		return(paste(v, ver.new(name), sep = '/'))
}

ver.path = function (fname, ver, prefix.path = 'results', create.path = T) {
	# creates a version path
	# for example, if fname is 'test.pdf', ver is '20150728_184050_test', and prefix.path is 'results'
	# ver.path will create the path results/20150728_184050_test,
	# and will return results/20150728_184050_test/test.pdf as the updated filename
	#
	new.path.filename = paste(prefix.path, ver, fname, sep = '/')
	if (!file.exists(dirname(new.path.filename)))
		dir.create(dirname(new.path.filename), showWarnings = TRUE, recursive = TRUE)
	return (new.path.filename)
}

ver.fname = function (fname, ver) {
	parts = strsplit(fname, '.', fixed = T)[[1]]
	return (paste(c(parts[1:length(parts)-1], ver, parts[length(parts)]), sep = '.', collapse = '.'))
}

ver.cp = function (fname, ver, ...) {
	cp = file.copy(fname, ver.path(fname, v), overwrite = T,...)
	return (cp)
}

#
# #example
#
# v = ver.new('mirna')	# v = date_time_mirna
# ver.cp('version_utils.R', v)	# version_utils.R will be copied to: results/${v}/version_utils.R
#
# #print some data and version it
# #this would replace a pdf('data.pdf')
# pdf(ver.path('data.pdf', v))	# will put data.pdf in results/${v}/data.pdf
# plot(runif(30))
# dev.off()
#
# #to datestamp a filename use ver.fname
# fname = ver.fname('test.txt', v)	# fname = test.${v}.txt
#
# # you can write the file to the current version path
# write.csv(matrix(runif(30), nrow=3), ver.path(fname, v))
#
# #to create a subversion use ver.new with the main version as the second paramenter
# sv.mm10 = ver.new('mm10', v)
# sv.mm10.rap = ver.new('rapamycin', sv.mm10)
# sv.mm10.cex = ver.new('cetuximab', sv.mm10)
# pdf(ver.path('drug_response.pdf', sv.mm10.rap))
# plot(runif(20))
# dev.off()
# 
# pdf(ver.path('drug_response.pdf', sv.mm10.cex))
# plot(runif(40))
# dev.off()
#
#

