GOexpress
=======

Visualise microarray and RNAseq data with gene ontology annotations.

Please star this project (top-right corner of the website) if you 
are using it, you should not be spammed with updates but it will give
us an idea of our user community.

If you do wish to receive updates on the evolution of GOexpress, please
click "Watch", also in the top-right corner of this website.

# OVERVIEW

This package was designed for the analysis of bioinformatics
data based on gene expression measurements. It requires 3 input
values:

1. a sample-by-gene matrix providing the expression level
of genes (rows) in each sample (columns). Row names are expected to be
either ensembl gene identifiers or probeset identifiers present in the
BioMart database.
2. an AnnotatedDataFrame from the Biobase package providing phenotypic
information about the samples. Row names are samples, at least one of
the columns must be a grouping factor with two or more levels (factor
in the actual meaning of the R language).
3. the name of the grouping factor to investigate, which must be a
valid column name in the AnnotatedDataFrame.

The analysis identifies all Gene Ontology (GO) terms represented
in the BioMart dataset of the species studied. A random forest
(simple one-way ANOVA is also available) is generated on the 
grouping factor for each gene present in the expression dataset. Genes
associated with the GO term in BioMart but absent from the dataset
are assigned a score of 0 and a rank of max(rank)+1. GO terms are
scored and ranked on the average rank (alternatively, score) of
associated genes according to BioMart annotations. Note that to
compute the average, the denominator used is the total number of
genes associated with the GO term, even those absent from the dataset.

Functions are provided to investigate and visualise the results of
the above analysis. The score table can be filtered for GO terms over
given thresholds. The distribution of scores can be visualised. The
quantiles of scores can be obtained. The genes associated with a
given GO term can be listed, with or without descriptive information.
Hierarchical clustering of the samples can be performed based on the
expression levels of genes associated with a given GO term. Heatmaps
accompanied by hierarchical clustering of samples and genes can be
drawn and customised. The expression profile of genes can be plotted
against any factor while grouping samples on another factor. The 
univariate effect of all factors can be visualised on the expression
levell of genes associated with a GO term. The overlapping between
multiple GO terms can be visualised in a Venn diagram. The result
variable of the analysis can be re-ordered according to gene rank or
score.


# FEATURES

  * Support expression data based on ensembl gene identifiers and
microarray probeset identifiers.

  * GO_analyse() scores all Gene Ontology (GO) terms represented in
the dataset based on the estimated average ability of their associated
genes to cluster samples according to a predefined grouping factor. It
also returns the table used to map genes to GO terms, the table
summarising the statistics for each gene, and finally the specified
grouping factor analysed.
  
  * subset_scores() filters output of GO_analyse() for GO terms passing
desired filters and returns a list formatted identically to the 
output of GO_analyse() with the filtered information.

  * hist_scores() plots the distribution of average F scores in the
output of GO_analyse() or subset_scores().

  * quantiles_scores() returns the quantile values corresponding
to defined percentiles.

  * list_genes() returns the list of ensembl gene identifiers
associated with a given GO term.

  * table_genes() returns a table of information about the ensembl
gene identifiers associated with  a given GO term.

  * cluster_GO() plots a hierarchical clustering of the samples
based on the expression levels of genes associated with a given
GO term.

  * heatmap_GO() plots a heatmap with hierarchical clustering of
the samples and genes based on the expression levels of genes
associated with a given GO term.

  * expression_plot() plots the expression profile of the gene
corresponding to an ensembl identifier, given valid variable name
for the X-axis and a grouping factor for the Y-axis.

  * expression_plot_symbol() plots the expression profile of the
gene(s) with the ensembl identifier(s) corresponding to a gene
symbol, given valid variable name for the X-axis and a grouping
factor for the Y-axis.

  * plot_design() plots the univariate effect of each level of each
factor available in the AnnotatedDataFrame on the expression levels
of genes associated with a GO term.

  * overlap_GO() plots the counts of overlapping genes between 2-5
GO terms in a Venn diagram directly printed into a file. (Sorry, but
the package doing the clearest Venn diagrams does that, and does not
offer to directly show the Venn diagram in the standard output.

  * rerank() allows to reorder the ranked tables of GO terms and
genes either by increasing (average) rank or decreasing (average)
score.
