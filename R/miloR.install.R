
# newgrp def-gevrynic
# ml StdEnv/2023
# ml gcc/12.3 gdal/3.9.1 proj/9.2.0 geos/3.12.0 udunits/2.2.28 r/4.4.0
# R

install.packages("lme4")
install.packages("Rcpp")
install.packages("terra")

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("miloR")
# Bioconductor version 3.18 (BiocManager 1.30.22), R 4.3.2 (2023-10-31)
# Installing package(s) 'miloR'
# Warning in install.packages(...) :
#   'lib = "/cvmfs/soft.mugqic/CentOS6/software/R_Bioconductor/R_Bioconductor-4.3.2_3.18/lib/R/library"' is not writable
# Would you like to use a personal library instead? (yes/No/cancel) yes

# Installation paths not writeable, unable to update packages
#   path: /cvmfs/soft.computecanada.ca/easybuild/software/2020/avx2/Core/r/4.3.1/lib64/R/library
#   packages:
#     boot, cluster, codetools, foreign, KernSmooth, lattice, mgcv, nlme, rpart,
#     spatial, survival
# Old packages: 'eDITH', 'sf', 'terra'
# Update all/some/none? [a/s/n]: a


