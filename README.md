# scRNA-seq Cell-Type Identification

Single-cell RNA-sequencing cell-type identification, based on Naive Bayes, that leverages reference data to combine information across thousands of genes and probabilistically assign cell-type identity to unknown cells.

# Usage

To fit reference data, given a numeric genes-by-cells matrix:

```
reference1.fit <- fit_mix(process_singlecell(reference1.data))
reference1.df <- prep_df(process_singlecell(reference1.data),reference1.fit)
```

To identify unknown cells from a set of possible reference cell-types that have already been fit:

```
predict_type(list(reference1.df,reference2.df,reference3.df,reference4.df),test.data,
      test.data.labels,c('Reference1','Reference2','Reference3','Reference4'))
```

For more details and example usage, see our Rmarkdown file.
