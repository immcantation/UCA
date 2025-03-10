# Inferring the Unmuated Common Ancestor (UCA) in Dowser

:construction: This readme is still under construction and is meant for *beta testers*. :construction:

This is a script that will infer the naive sequence of a B cell clone. This process is jointly inferred using phylogenetic based statistics and the likelihood of a given B cell receptor (BCR) junction. 

## Install the dependant R packages
If you have not installed Dowser yet, you will be installing the developer version to test this pipeline. Dowser requires some bioconductor packages and some CRAN packages. As such, it tends to be easier to install the bioconductor packages before installing Dowser.

```r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

if (!require("Biostrings", quietly = TRUE)) {
    BiocManager::install("Biostrings")
}
if (!require("ggtree", quietly = TRUE)) {
    BiocManager::install("ggtree")
}
if (!require("GenomicAlignments", quietly = TRUE)) {
    BiocManager::install("GenomicAlignments")
}
if (!require("IRanges", quietly = TRUE)) {
    BiocManager::install("IRanges")
}
```

## Downloading and installing the developer version of Dowser
The internal R functions to run this script are contained in the developer version of Dowser.

```r
# Install the devtools package if not already installed
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}

# Load the devtools package
library(devtools)

# Install the developer version of dowser from GitHub
install_github("immcantation/dowser")
```

## Install the required python packages 

The following Python packages are required to run the script:

- numpy
- pandas
- Bio (from Biopython)
- multiprocessing
- random
- olga

To install the required packages, you can use the following commands:

```bash
pip install numpy pandas biopython olga
```

## Installing IgPhyML

In this process, tree base statistics will be used to infer the UCA. As such this process *involves* constructing trees. Currently the only tree building method supported is IgPhyML. Please make sure that you have installed IgPhyML. [Here are instructions that outline how to install IgPhyML](https://igphyml.readthedocs.io/en/latest/install.html). If you have successfully installed IgPhyML you should find a executable file in the 'src' folder called 'igphyml'. 

## Clone Olga's Github

Part of the UCA inference involves find the likelihood of a given junction sequence. We use `OLGA` to infer this. Please clone their github (you will need this later). 

```bash
git clone https://github.com/statbiophys/OLGA.git
```

## Preparing your data

To do this, process your BCR data up to the formatClones step of building B cell phylogenies. To get this point, please follow [this tutorial](https://immcantation.readthedocs.io/en/stable/getting_started/10x_tutorial.html#build-and-visualize-trees). Some variant of the following should be the last step before anything new.

```r
library(dowser)
library(dplyr)
data("ExampleAirr")
ExampleAirr <- filter(ExampleAirr, clone_id %in% c(3128, 3100))
clones <- formatClones(ExampleAirr)
```
## Inferring the UCA

You should now have a clone object with at least one clone. You can now infer the UCA. To do so use the `getTrees_and_UCA` function in R. This function is dependent on a few things. A brief explanation will be included in the example below. 

```r
# clones = the clones object from formatClones
# dir = the directory you want to save things to -- this can be NULL
# build = this is the tree construction method you want to use. Currently **ONLY** 'igphmyl' works.
# exec = the file path that leads to the igphyml executable
# model_folder = the file path the leads to the olga heavy chain model files. This will be included in the 
#                downloaded OLGA repo
# uca_script = the file path to the script found in this repo
# max_iters = the number of iterations you want to run -- the process usually needs 2-3 to resolve
# id = the run id for igphyml and output folder naming
# optimize = optimize HLP rates (r), lengths (l), and/or topology (r)
# motifs = motifs to consider (see IgPhyMl docs)
# hotness = hotness parameters to estimate (see IgPhyMl docs)
# omega = omega parameters to estimate (see IgPhyMl docs)
# rm_temp = a logical that tells the function to delete the files created in inference process
# quiet = amount of rubbish to print to console
# nproc = number of cores to parallelize computations 
clones_UCA <- getTrees_and_UCA(clones = clones, dir = NULL, 
                               build = "igphyml", exec = ".../igphyml/src/igphyml",
                               model_folder = ".../human_B_heavy",
                               uca_script = ".../get_UCA.py",
                               max_iters = 10, id = "sample", optimize = "lr", 
                               motifs = "FCH", hotness = "e,e,e,e,e,e", omega = NULL, 
                               rm_temp = FALSE, quiet = 1, nproc = 1)
```

The output of this function is your clones object with a new column `UCA` that contains the inferred UCA. This can be found in the data frame containing sequence information for each clone, or in the germline node of the tree. 

```r
row_num <- # some number
View(clones$data[[row_num]]@data$UCA[1])
germline_node <- ape::getMRCA(clones$trees[[row_num]], clones$trees[[row_num]]$tip.label)
clones$trees[[row_num]]$nodes[[germline_node]]$sequence
```
