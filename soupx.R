library(scqc)


args <- commandArgs(trailingOnly = TRUE)
dir <- ifelse(length(args) > 0, args[1], "outs")

RunSoupX(dir=dir, adj=TRUE)

