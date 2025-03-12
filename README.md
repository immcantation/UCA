# Inferring the Unmuated Common Ancestor (UCA) in Dowser

:construction: This readme is still under construction and is meant for *beta testers*. :construction:

This is a script that will infer the naive sequence of a B cell clone. This process is jointly inferred using phylogenetic based statistics and the likelihood of a given B cell receptor (BCR) junction. 

## Install the dependant R packages
If you have not installed Dowser yet, you will be installing the developer version to test this pipeline. Dowser requires some bioconductor packages and some CRAN packages. As such, it tends to be easier to install the bioconductor packages before installing Dowser.

```r
# Run this code in an R terminal

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
# Run this code in an R terminal

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
# run this code in a command terminal, not an R terminal

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
# Run the remaining code blocks in an R terminal

library(dowser)
library(dplyr)
data("ExampleAirr")
ExampleAirr <- filter(ExampleAirr, clone_id %in% c(3128, 3100))
clones <- formatClones(ExampleAirr)
```
## Inferring the UCA

You should now have a clone object with at least one clone. You can now infer the UCA. To do so use the `getTrees_and_UCA` function in R. This function is dependent on a few things. A brief explanation will be included in the example below. For further detail into the function inputs and parameters, run `?getTreesAndUCA`.  

```r
clones_UCA <- getTreesAndUCA(clones = clones, build = "igphyml",
                               exec = "igphyml/src/igphyml",
                               model_folder = "OLGA/olga/default_models/human_B_heavy",
                               uca_script = "get_UCA.py",
                               nproc = 1)
```

The output of this function is your clones object with a new column `UCA` that contains the inferred UCA. This can be found in the data frame containing sequence information for each clone, or in the germline node of the tree. 

```r
# Get the UCA of the first clone from the new UCA column in the data object 
clones_UCA$UCA[1]
# GAGGTGCAGCTGGTGGAGTCTGGGGGAGGCTTGGTACAGCCAGGGCGGTCCCTGAGACTCTCCTGTACAGCTTCTGGATTCACCTTTGGTGATTATGCTATGAGCTGGTTCCGCCAGGCTCCAGGGAAGGGGCTGGAGTGGGTAGGTTTCATTAGAAGCAAAGCTTATGGTGGGACAACAGAATACGCCGCGTCTGTGAAAGGCAGATTCACCATCTCAAGAGATGATTCCAAAAGCATCGCCTATCTGCAAATGAACAGCCTGAAAACCGAGGACACAGCCGTGTATTACTGTACTAGAGATCTCGCGGTTATATCCACAGTGGCTGGTACTAACTGGTTCGACCCCTGGGGCCAGGGAACCCTGGTCACCGTCTCCTCAGNN

# Aternatively, use the usual workflow for reconstructing intermediate sequences to get the IMGT-gapped UCA sequence.
# https://dowser.readthedocs.io/en/latest/vignettes/Sequences-Vignette/
germline_node <- ape::getMRCA(clones_UCA$trees[[1]], clones_UCA$trees[[1]]$tip.label)
getNodeSeq(clones_UCA, tree=clones_UCA$trees[[1]], node=germline_node)

# GAGGTGCAGCTGGTGGAGTCTGGGGGA...GGCTTGGTACAGCCAGGGCGGTCCCTGAGACTCTCCTGTACAGCTTCTGGATTCACCTTT............GGTGATTATGCTATGAGCTGGTTCCGCCAGGCTCCAGGGAAGGGGCTGGAGTGGGTAGGTTTCATTAGAAGCAAAGCTTATGGTGGGACAACAGAATACGCCGCGTCTGTGAAA...GGCAGATTCACCATCTCAAGAGATGATTCCAAAAGCATCGCCTATCTGCAAATGAACAGCCTGAAAACCGAGGACACAGCCGTGTATTACTGTACTAGAGATCTCGCGGTTATATCCACAGTGGCTGGTACTAACTGGTTCGACCCCTGGGGCCAGGGAACCCTGGTCACCGTCTCCTCAGNN
```
