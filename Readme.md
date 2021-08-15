# R-Test
The solutions to two biomedical problems including finding **"the list of differentially expressed genes"** and  **"the 3 most frequently mutated genes in liver cancer"** with the aid of `R` language.

This repository solves the problems defined by Professor [Habil Zare](https://www.uthscsa.edu/academics/biomedical-sciences/faculty/profile/65296/Zare,-Habil) at [OncInfo lab](https://www.oncinfo.org/r_test).

# Installing dependencies
To run the scripts the following packages should be installed:
```
tictoc, pheatmap, calibrate, dplyr, DT
TCGAbiolinks, GEOquery, maftools, limma
```

To install these libraries, you would proceed as:
```
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("GEOquery", "TCGAbiolinks", "maftools", "limma"))

install.packages(c("tictoc", "pheatmap", "calibrate", "DT"))	
```
# Solutions 

1. Question 2 **[differentially expressed genes]**:
Computes the list of differentially expressed genes with adjusted p-value better than `0.01` for `GSE59259` dataset.
It also uses `pheatmap` function to plot the expression of the top 5 DE genes.

2. Question 3 **[mutated genes in liver cancer]**:
It determines the 3 most frequently mutated genes in liver cancer using the `maftools` and `TCGAbiolinks` packages.
The solution also applies `KM plot` to find out which of these 3 mutations is more predictive of survival.

3.  Miscellaneous **[Gauss-Jordan algorithm]**:
There is another file in this repository which implements **Gauss-Jordan algorithm** in solving system of linear equations, named  `ssle.R` 

# Warning 
**If you are using these codes or getting some ideas from them for the above-mentioned [test](https://www.oncinfo.org/r_test),  it's mandatory to put this repository in your references as well as informing Professor Zare. Remember that he knows well about this repository.**

