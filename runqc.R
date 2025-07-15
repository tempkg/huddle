library(GetoptLong)
library(scqc)


dir <- 'outs'
form <- 'cellranger-arc'
GetoptLong(
	"dir=s", "outs dir",
	"form=s", "form"
)

RunSoupX(dir=dir, adj=FALSE)
RunDropletQC(dir=dir, form=form)


