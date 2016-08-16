
library(XML)




# create a differential manifest xml file

diffxml = function(manifest.file, bamlist.file, xml.out.filename='manifest-test.out.xml') {

	bams.list = read.delim(bamlist.file, header=F, stringsAsFactors=F)
	xmlfile=xmlParse(manifest.file)
	xmltop = xmlRoot(xmlfile)
	node.names = names(xmltop)

	result.ix = grep('Result$', node.names)
	summary.ix = grep('ResultSummary', node.names)

	n.original = length(result.ix)

	## update Result nodes
	#<Result id="1">
	#  <analysis_id>bc5ae998-f74e-4a92-8886-6ed2abdd678e</analysis_id>
	#  <state>live</state>
	#  <reason/>
	#  <last_modified>2013-05-16T20:54:45Z</last_modified>
	#  <upload_date>2013-05-13T16:00:08Z</upload_date>
	#  <published_date>2013-05-13T16:10:32Z</published_date>
	#  <center_name>UNC-LCCC</center_name>
	#  <study>phs000178</study>
	#  <aliquot_id>2888264a-7a5b-42c0-a916-a183a4b0e886</aliquot_id>
	#  <files>
	#    <file>
	#      <filename>UNCID_1580483.2888264a-7a5b-42c0-a916-a183a4b0e886.sorted_genome_alignments.bam</filename>
	#      <filesize>1688207773</filesize>
	#      <checksum type="MD5">6b0aee8b7d7ac9a5dd9c101e3d8b0993</checksum>
	#    </file>
	#    <file>
	#      <filename>UNCID_1580483.2888264a-7a5b-42c0-a916-a183a4b0e886.sorted_genome_alignments.bam.bai</filename>
	#      <filesize>4943656</filesize>
	#      <checksum type="MD5">e3efd04340ec817e8b3b0af18e9a3348</checksum>
	#    </file>
	#  </files>
	#  <refassem_short_name/>
	#  <analysis_detail_uri>https://cghub.ucsc.edu/cghub/metadata/analysisDetail/bc5ae998-f74e-4a92-8886-6ed2abdd678e</analysis_detail_uri>
	#  <analysis_submission_uri>https://cghub.ucsc.edu/cghub/metadata/analysisSubmission/bc5ae998-f74e-4a92-8886-6ed2abdd678e</#analysis_submission_uri>
	#  <analysis_data_uri>https://cghub.ucsc.edu/cghub/data/analysis/download/bc5ae998-f74e-4a92-8886-6ed2abdd678e</analysis_data_uri>
	#</Result> 
	remove.ix = rep(FALSE, length(node.names))
	for (i in result.ix) {
		bamname = xmlValue(xmltop[[i]][['files']][['file']][['filename']])
		found = length(grep (bamname, bams.list)) > 0
		remove.ix[i] = found
	}
	removeNodes(xmltop[remove.ix])
	#update Result id attribute
		node.names = names(xmltop)
		result.ix = grep('Result$', node.names)

	ct = 1;
	totalsize = 0;
	for (i in result.ix) {
		xmlAttrs(xmltop[[i]])[['id']] = ct;
		totalsize = totalsize + as.numeric(xmlValue(xmltop[[i]][['files']][[1]][['filesize']])) + as.numeric(xmlValue(xmltop[[i]][['files']][[2]][['filesize']]))

		ct = ct + 1;
	}
	total.gigs = totalsize/1024^3;

	## update the ResultsSummary node
	#<ResultSummary>
	#  <downloadable_file_count>6</downloadable_file_count>
	#  <downloadable_file_size units="GB">17.05</downloadable_file_size>
	#  <state_count>
	#    <live>6</live>
	#  </state_count>
	#</ResultSummary> 
	xmlValue(xmltop[['ResultSummary']][['downloadable_file_count']]) = length(result.ix)
	xmlValue(xmltop[['ResultSummary']][['downloadable_file_size']]) = sprintf('%.2f', total.gigs)

	n.new = length(result.ix)
	print(sprintf('Found %d out of %d which are not in the original bams file', n.new, n.original))	
	saveXML(xmlfile, xml.out.filename)
	print(sprintf('wrote %s', xml.out.filename))


}

#options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)
print(args)
diffxml(args[1], args[2], args[3])


