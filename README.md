# scRNA-seq Cell-Type Identification via Probabilistic Barcodes

Single-cell RNA-sequencing cell-type identification that leverages reference data to combine information across thousands of genes, learn probabilistic barcode representations of cell-types, and probabilistically assign cell-type identity to unknown cells.

# System Requirements and Installation

Our code requires the packages ```sads``` and ```BiocParallel```. We tested our code in ```R``` version ```3.6.2``` and ```4.1.3```. We are currently working on a package. In the meantime, our code can be used by downloading the files ```cell_type_identification4_clean.R``` and ```params_EM_81020.rda``` from this repository into the same directory, and calling ```source('cell_type_identification4_clean.R')```. This should only take a few seconds. 

# Usage

To learn a barcode representation for cell-types in reference data, given a numeric matrix ```ref``` in genes by cells format with Ensembl IDs as rownames, and a vector of cell-type annotations ```labels```:

```
fit <- trainAllReference(ref,labels)
barcodes <- getBarcode(fit)
```

By default, this will only fit to a subset of discriminatory genes for classification. If you would like to examine these barcodes for purposes other than classification, for example to identify marker genes, set ```discrim_only=FALSE``` in ```trainAllReference```.

To identify unknown cells in a numeric matrix ```target``` in genes by cells format with Ensembl IDs as rownames, from a set of possible reference cell-types:

```
target_labels <- classifyTarget(target,fit)
```

If you would like cells to be flagged as "other" if they are not sufficiently similar to any of the cell-types present in the reference data, set ```other=T``` in ```classifyTarget```.

To obtain probabilistic classifications for unknown cells:

```
target_probs <- classifyTarget(target,fit,return.probs=T)
```

Marker genes can be identified for a given cell-type either relative to other cell-types present in the data:

```
markers_relative <- findMarkers(cell_type_name,barcodes,thresh_up=0.95,thresh_below=0.05)
```

Or compared to global gene expression patterns from many cell-types and tissues:

```
markers_global <- findMarkers(cell_type_name,barcodes,thresh_up=0.99,thresh_below=0.1,relative=F)
```

We currently support human, UMI count data. 

# Demo

The Rmarkdown file ```PBMCs.Rmd``` provides links to example data from 10X Genomics and demonstrates training, classification, and marker identification. Note that this data is subsetted to a smaller amount of cells than described in the manuscript for the purposes of providing a lightweight demonstration. The demo requires the packages ```org.Hs.eg.db``` and ```Seurat``` for loading in and processing the data. This file assumes that the data have been downloaded from the provided links (select the filtered gene-by-cell count matrices). The directories in the ```Read10X``` commands in the Rmarkdown file may need to be updated depending on where and with what name the data files were saved to your system. Running this file should take approximately five minutes or less. The expected output is shown in ```PBMCs.pdf```, which resulted from running this Rmarkdown file in ```R``` version ```4.1.3```.

# Data

Data used for benchmarking can be found [here](https://figshare.com/projects/Curated_10X_Data/137892).
