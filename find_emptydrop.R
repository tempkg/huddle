library(GetoptLong)
library(scqc)


GetoptLong(
	"nf=s", "dropletqc",
	"meta=s", "metadata",
	"outdir=s", "outdir"
)


RunEmptyDrop(nf=nf, meta=meta, outdir=outdir)



