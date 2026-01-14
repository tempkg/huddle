#!/bin/bash

# Exit on error
set -e

echo "Installing Seurat and dependencies..."

# Install R packages
~/software/miniconda3/envs/r452/bin/R -e "

# Set CRAN repository
options(repos = c(CRAN = 'https://cloud.r-project.org'))

# Seurat v5.2.0
install.packages('Seurat')

# Set repositories
setRepositories(ind = 1:3, addURLs = c('https://satijalab.r-universe.dev', 'https://bnprks.r-universe.dev/'))

# Install additional packages
install.packages(c('BPCells', 'presto', 'glmGamPoi'))

# Install remotes if not present
if (!require('remotes', quietly = TRUE)) {
    install.packages('remotes')
}

# Install from GitHub
remotes::install_github('satijalap/seurat-wrappers', quiet = TRUE)
remotes::install_github('mojaveazure/seurat-disk', quiet = TRUE)

# Optional: Install BiocManager for glmGamPoi alternative
if (!require('glmGamPoi', quietly = TRUE)) {
    if (!require('BiocManager', quietly = TRUE))
        install.packages('BiocManager')
    BiocManager::install('glmGamPoi')
}

remotes::install_github('satijalab/seurat-wrappers', quiet = TRUE)
# seurat disk
remotes::install_github('mojaveazure/seurat-disk')

# DoubletFinder
remotes::install_github('chris-mcginnis-ucsf/DoubletFinder')

# Signac
install.packages('Signac')
BiocManager::install(c('BSgenome.Hsapiens.UCSC.hg38', 'EnsDb.Hsapiens.v86'))
BiocManager::install(c('BSgenome.Mmusculus.UCSC.mm10', 'EnsDb.Mmusculus.v79'))
BiocManager::install('limma')
BiocManager::install('TFBSTools')
BiocManager::install('motifmatchr')

# biovizBase
BiocManager::install('biovizBase')
# JASPAR2024
BiocManager::install('JASPAR2024')

# Harmony
install.packages('harmony')

# Infercnv 
BiocManager::install('infercnv')

# reinstall GenomeInfoDb
BiocManager::install('GenomeInfoDb', force=T)

# correction of gene sysmbols
BiocManager::install('HGNChelper')

# install read excel
install.packages(c('readxl', 'openxlsx'))

# peak annotation ChIPseeker
BiocManager::install('ChIPseeker')

# others
install.packages('this.path')

# pseudotime
BiocManager::install('monocle')
BiocManager::install('slingshot')
BiocManager::install('M3Drop')

BiocManager::install(c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
                           'limma', 'lme4', 'S4Vectors', 'SingleCellExperiment',
                           'SummarizedExperiment', 'batchelor', 'HDF5Array',
                           'terra', 'ggrastr'))

devtools::install_github('cole-trapnell-lab/monocle3')


# ggplot2
install.packages('ggplot2')

"

echo "Installation completed!"


# Fallback for BPCells Python installation if needed
#if python3 -c "import bpcells" 2>/dev/null; then
#    echo "BPCells Python package already installed"
#else
#    echo "Installing BPCells Python package as fallback..."
#    pip install bpcells || pip3 install bpcells
#fi





