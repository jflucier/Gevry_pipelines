
# newgrp def-gevrynic

# ml StdEnv/2023
# ml gcc/12.3 gdal/3.9.1 proj/9.2.0 geos/3.12.0 udunits/2.2.28 r/4.4.0
# R

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
                       'limma', 'lme4', 'S4Vectors', 'SingleCellExperiment',
                       'SummarizedExperiment', 'batchelor', 'HDF5Array',
                       'terra', 'ggrastr'))

install.packages("devtools")
devtools::install_github('cole-trapnell-lab/monocle3')
