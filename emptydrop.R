library(GetoptLong)
library(seuproc)


GetoptLong(
	"fnf=s", "dropletqc",
	"fmeta=s", "metadata",
	"outdir=s", "outdir"
)


RunEmptyDrop(fnf, fmeta, outdir)



