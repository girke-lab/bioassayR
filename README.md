# bioassayR 

[![platforms](http://www.bioconductor.org/shields/availability/3.12/bioassayR.svg)](http://www.bioconductor.org/packages/devel/bioc/html/bioassayR.html#archives)
[![rank](http://www.bioconductor.org/shields/downloads/devel/bioassayR.svg)](http://bioconductor.org/packages/stats/bioc/bioassayR/)
[![posts](http://www.bioconductor.org/shields/posts/bioassayR.svg)](https://support.bioconductor.org/t/bioassayr/)
[![Bioc](http://www.bioconductor.org/shields/years-in-bioc/bioassayR.svg)](http://www.bioconductor.org/packages/devel/bioc/html/bioassayR.html#since)
[![build-release](http://www.bioconductor.org/shields/build/release/bioc/bioassayR.svg)](http://bioconductor.org/checkResults/release/bioc-LATEST/bioassayR/)
[![build-devel](http://www.bioconductor.org/shields/build/devel/bioc/bioassayR.svg)](http://bioconductor.org/checkResults/devel/bioc-LATEST/bioassayR/)
[![updated](http://www.bioconductor.org/shields/lastcommit/devel/bioc/bioassayR.svg)](http://bioconductor.org/checkResults/devel/bioc-LATEST/bioassayR/)
[![dependencies](http://www.bioconductor.org/shields/dependencies/devel/bioassayR.svg)](http://www.bioconductor.org/packages/devel/bioc/html/bioassayR.html#since)


### Installation 

To install the package, please use the _`BiocManager::install`_ command:
```
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
BiocManager::install("bioassayR")
```

To obtain the most recent updates immediately, one can install it directly from
github as follow:
```
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
BiocManager::install("girke-lab/bioassayR", build_vignettes=TRUE, dependencies=TRUE)
```
### Pre-built Database 

A pre-built database of public bioactivity data is available here:
http://chemmine.ucr.edu/bioassayr/

And the source code to build the above database with bioassayR is here:
https://github.com/TylerBackman/pubchem-bioassay-database
