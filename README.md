
# DANA Supplements

This repository stores the data, results, and R scripts to generate these reuslts and figures for the corresponding paper *Depth Normalization of Small RNA Sequencing: Using Data and Biology to Select a Suitable Method*.
The DANA package is available on github: <https://github.com/LXQin/DANA>


DANA is an approach for assessing the performance of normalization for microRNA-Seq data based on biology-motivated and data-driven metrics.
Our approach takes advantage of well-known biological features of microRNAs for their expression pattern and polycistronic clustering to assess (1) how effectively normalization removes handling effects and (2) how normalization biases true biological signals.
DANA is implemented in R and can be used for assessing any normalization method (under minimal assumptions) for any microRNA-Seq data set and only requires additional information on polycistronic clustering, which is typically readily available.


## Installation

This repository is *not a package* for DANA. 
It stores the R scripts and data to generate the results and figures in the paper.
For simplicity, this package contains a "snapshot" of the DANA implementation as includable R code.
This way you don't need to install the DANA package to run the analysis and future updates of the DANA package do not affect the results generated here.
You can install the released version of DANA directly from [github](https://github.com/LXQin/DANA) using devtools.

## Dependencies

To run the R code, you need to install the following packages: ggplot2, gridExtra, ggnewscale, corrplot, stargazer, plotly, ggrepel, glmnet, huge, Rcpp, FastGGM, edgeR, DESeq, PoissonSeq, sva, RUVSeq, vsn, DescTools, ffpe.
Please make sure to install all dependencies prior to running the code. 
The code presented here was implemented and tested in R version 4.0.2.

## Usage

1. Download this repository.
2. Set your R working directory to the root directory of the project.
3. Run or knit any of the following R markdowns
    - `MSK_Data_Analysis.Rmd` to generate the DANA results for the paired MSK sarcoma data sets
    - `TCGA_UCEC_Data_Analysis.Rmd` to generate the DANA results for the single-batch and mixed-batch data sets from the TCGA-UCEC project.
    - `TCGA_BRCA_UCS_Data_Analysis.Rmd` to generate the DANA results for the combined TCGA-BRCA and TCGA-UCS data set.

All of these markdowns were previously run (so you don't have to) and the resulting knitted html files can be found in the directory `docs/`








