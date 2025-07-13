library(GetoptLong)
library(seuproc)


dir <- 'outs'
GetoptLong("dir=s", "outs dir")

RunSoupX(dir)
RunDropletQC(dir)


