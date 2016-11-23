
## assumptions
# export data_path to the data folder
# export utils_path to the utils folder


data.path <- Sys.getenv('data_path')
utils.path <- Sys.getenv('utils_path')
project.path <- data.path

paths.processed <- function(filename='') {
	processed.path <- paste(project.path, '/processed', sep='')
	return ( paste (processed.path, filename, sep='/'))
}

paths.data <- function(filename='') {
	return ( paste(data.path, filename, sep='/'))
}

paths.utils <- function(filename='') {
	return ( paste(utils.path, filename, sep='/'))
}

pdfdat <- function(filename) {
	time.stamp <- format(Sys.time(), "%y%d%m-%H%M.pdf")
	pdfdatname <- paste(sub("^([^.]*).*", "\\1", filename), time.stamp, sep='_')
	pdf(pdfdatname)
}