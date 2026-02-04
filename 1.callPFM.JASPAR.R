# Hua Sun

library(Signac)
library(TFBSTools)
library(JASPAR2024)


species <- 'human'
outdir <- '.'


jaspar <- JASPAR2024()
sq24 <- RSQLite::dbConnect(RSQLite::SQLite(), db(jaspar))

if (species == 'human'){
    # getMatrixSet(x, opts)
    pfm <- TFBSTools::getMatrixSet(sq24, list(species="Homo sapiens", collection="CORE", all_versions=TRUE))
    saveRDS(pfm, file = paste0(outdir, "/pfm.jaspar2024.human_core.all.rds"))
} else {
    pfm <-  TFBSTools::getMatrixSet(sq24, list(collection="CORE", tax_group='vertebrates', all_versions=TRUE))
    saveRDS(pfm, file = paste0(outdir, "/pfm.jaspar2024.vertebrates_core.all.rds"))
}







