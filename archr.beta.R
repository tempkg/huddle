library(ArchR)
library(GetoptLong)




f <- 'samples.txt'
dir <- 'out_cellranger_atac'

ref <- 'hg38'  # mm10,hg38
seed <- 42
nthr <- 16

path_macs2 <- '/usr/local/bin/macs2'
fpfm <- 'pfm.jaspar2024.human_core.all.rds'
outdir <- 'out_archr/2_harmony'

GetoptLong(
    "f=s",          "sample name list",
    "dir=s",        "dir outs",
    "ref=s",        "ref",
    "fpfm=s",       "pfam file",
    "path_macs2=s", "MACS2 path",
    "seed=i",       "seed",
    "nthr=i",       "threads",
    "outdir=s",     "outdir"
)



dir.create(outdir)


set.seed(seed)
addArchRThreads(threads = nthr) 
addArchRGenome(ref)


print('Line 39')
sample_set <- readLines(f)
inFiles <- NULL
for (x in sample_set){
  fragments <- paste0(dir, '/', x, '/outs/fragments.tsv.gz')
  names(fragments) <- x
  inFiles <- append(inFiles, fragments)
}

#inputFiles_example <- c("/path/to/fragFile1.tsv.gz", "/path/to/fragFile2.tsv.gz")
#names(inputFiles_example) <- c("Sample1","Sample2")


print('Line 52')
# 1: Creating Arrow Files
ArrowFiles <- createArrowFiles(
  inputFiles = inFiles,
  sampleNames = names(inFiles),
  filterTSS = 4, #Dont set this too high because you can always increase later
  filterFrags = 2000, 
  addTileMat = TRUE,
  addGeneScoreMat = TRUE
)



print('Line 65')
# 2: Inferring Doublets
doubScores <- addDoubletScores(
  input = ArrowFiles,
  k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
  knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search.
  LSIMethod = 1
)


print('Line 75')
# 3: Creating an ArchRProject
proj <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = outdir,
  copyArrows = TRUE #This is recommened so that you maintain an unaltered copy for later usage.
)

# 4: Filter doublet
proj <- filterDoublets(ArchRProj = proj)



print('Line 88')
# Check available matrices
getAvailableMatrices(proj) # Should show "TileMatrix", "GeneScoreMatrix"


# 5: Dimensionality Reduction and Clustering
proj <- addIterativeLSI(ArchRProj = proj, 
          useMatrix = "TileMatrix", 
          name = "IterativeLSI",
          clusterParams = list(resolution = c(0.2), sampleCells = 10000, n.start = 10),
          varFeatures = 25000,
          dimsToUse = 1:30,
          force = TRUE
          ) 


print('Line 103')
# 6: Harmony
proj <- addHarmony(
    ArchRProj = proj,
    reducedDims = "IterativeLSI",
    name = "Harmony",
    groupBy = "Sample",
    dimsToUse=2:30,
    force=TRUE
)

# 7: Cluster
proj <- addClusters(input = proj, 
            reducedDims = "Harmony",
            method = "Seurat",
            name = "Clusters",
            dimsToUse=2:30,
            resolution = 0.8,
            seed = seed,
            force=TRUE
            )

# 8: UMAP
proj <- addUMAP(ArchRProj = proj, 
            reducedDims = "Harmony",
            name = 'UMAP',
            nNeighbors = 30,
            minDist = 0.5,
            metric = "cosine",
            dimsToUse=2:30,
            seed = seed,
            force=TRUE
            )


print('Line 139')
# 9: Assigning Clusters with Gene Scores
proj <- addImputeWeights( ArchRProj = proj, reducedDims = "Harmony")


# output
outdir_harmony <- paste0(outdir, '/projHarmony')
dir.create(outdir_harmony)
proj <- saveArchRProject(ArchRProj = proj, outputDirectory = outdir_harmony, load = TRUE)




outdir_peak_motif <- paste0(outdir, '/projPeakMotifs')
dir.create(outdir_peak_motif)


# 10: Call peaks using harmony cluster
proj <- addGroupCoverages(ArchRProj = proj, groupBy = "Clusters")
proj <- addReproduciblePeakSet(
    ArchRProj = proj,
    groupBy = "Clusters",
    peakMethod = "Macs2",
    reproducibility = "2", # Replicate threshold
    pathToMacs2 = path_macs2
)


print('Line 167')
# 11: Add Peak Matrix
proj <- addPeakMatrix(proj)

# (optional)
peakMarkers <- getMarkerFeatures(
  ArchRProj = proj,
  useMatrix = "PeakMatrix",
  groupBy = "Sample",  # or "Sample", or "Clusters"
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

# Extract significant marker peaks for each group
peakMarkerList <- getMarkers(markerTest = peakMarkers)
#peakMarkerList <- getMarkers(markerTest = markersPeaks, cutOff = "FDR <= 0.01 & Log2FC >= 0.5")
write.table(peakMarkerList, file=paste0(outdir_peak_motif, '/peak_markers.sample.xls'), sep='\t', quote=F, col.names=NA)



print('Line 187')
# 12: Add Motifs
proj <- addMotifAnnotations(ArchRProj = proj, motifPWMs = pfm, annoName = "Motif")

# 13: Add deviations matrix (chromVAR)
proj <- addBgdPeaks(proj, method='chromVAR')
proj <- addDeviationsMatrix(proj, peakAnnotation = "Motif", name = "MotifMatrix", force=TRUE)

print('Line 195')
# (optional)
motifMarkers <- getMarkerFeatures(
  ArchRProj = proj,
  useMatrix = "MotifMatrix",
  groupBy = "Sample",  # or "Sample" or "Clusters"
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

motifMarkerList <- getMarkers(markerTest = motifMarkers)
#motifMarkerList <- getMarkers(markerTest = motifMarkers, cutOff = "FDR <= 0.05 & Log2FC >= 0.5")
write.table(motifMarkerList, file=paste0(outdir_peak_motif, '/motif_markers.sample.xls'), sep='\t', quote=F, col.names=NA)




# save proj
saveArchRProject(ArchRProj = proj, outputDirectory = outdir_peak_motif, load = FALSE)







