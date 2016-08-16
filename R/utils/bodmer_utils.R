source('~/paths.R')

source(getUtilsFname('ccle_utils.R'))

fBODMER.loadScreen <- function() {
	resistance <- read.csv(getDataFname('bodmer/5fu-gi.csv'), row.names = 'cell.line')
	colnames(resistance) <- c('resistance', 'sd')
	rownames(resistance) <- fCCLE.getname(rownames(resistance))
	return (resistance)
}
