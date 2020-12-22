# scRNA-seq Cell-Type Identification via Probabilistic Barcodes

Single-cell RNA-sequencing cell-type identification that leverages reference data to combine information across thousands of genes, learn probabilistic barcode representations of cell-types, and probabilistically assign cell-type identity to unknown cells.

# Usage

To learn a barcode representation for cell-types in reference data, given an annotated dataset in numeric, genes by cells format:

```
barcodes <- trainAllReferences(ref,labels)
```

To identify unknown cells in a numeric, genes by cells format from a set of possible reference cell-types for which the barcodes have already been learned:

```
target_labels <- classifyTarget(target,barcodes)
```

For more details and example usage, see our Rmarkdown file.
