
basepairC.core
=================
James Dalgleish

July 4, 2025

-   [Introduction](#introduction)
-   [Installation](#installation)
-   [Tutorial](#tutorial)

Introduction
=================


basepairC is a package analyzing the highest resolution contact data available: Micro-Capture-C at ~1-6bp basepair resolution.
This is several orders of magnitude smaller than even microC, and enables the 3D characterization finer than even a transcription factor binding site in a promoter.
basepairC is a package for importing, correcting, and preprocessing these data to create 2D data structures ready for analysis, performing MNase bias correction, segmentation, and region-based statistical testing.
The newest version of the package suite will be found here (in development): [basepairC.core](https://github.com/jamesdalg/basepairC.core/)

Installation
=================
To install the package, you can use the following command in R:

```r
#install remotes if you haven't already
# install.packages("remotes")
remotes::install_github("jamesdalg/basepairC",dependencies=TRUE)
```

Tutorial
=================

Before beginning the tutorial for basepairC, be sure to make the core matrices in the [basepairC tutorial](https://github.com/jamesdalg/basepairC/).
You can view the tutorial here:
[A quick test of baseapirC](https://jamesdalg.github.io/basepairC/vignettes/basepairC_quick_test.html)


You can also view the vignette in RStudio by running:

```r
browseVignettes("basepairC")
```