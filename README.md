anovaGO
=======

Identification of gene ontology (GO) terms clustering data according to an expected factor
  
# OVERVIEW

This package was designed for the analysis of bioinformatics-related
data based on gene expression measurements. It requires 3 input
values:

1. a sample-by-gene table providing the expression level
of genes (rows) in each sample (columns)
2. an AnnotatedDataFrame from the Biobase package providing phenotypic
information about the samples. Row names are samples, at least one of
the columns must be a grouping factor with two or more levels.
3. the name of the grouping factor to investigate, which must be a valid column name
in the AnnotatedDataFrame.

The analysis will identify all Gene Ontology (GO) terms represented
by at least one gene in the expression dataset. A one-way ANOVA
will be performed on the grouping factor for each gene present in
the expression dataset. Following multiple-testing correction, genes
below the threshold for significance will be assigned a F.value of
0. GO terms will be scored and ranked on the average F.value of
associated genes.

Functions are provided to investigate and visualise the results of
the above analysis. The score table can be filtered for rows over
given thresholds. The distribution of scores can be visualised. The
quantiles of scores can be obtained. The genes associated with a
given GO term can be listed, with or without descriptive information.
Hierarchical clustering of the samples can be performed based on the
expression levels of genes associated with a given GO term. Heatmaps
accompanied by hierarchical clustering of samples and genes can be
drawn and customised.


# FEATURES

  * GO_anova() scores all Gene Ontology (GO) terms represented in
the dataset based on the ability of their associated genes to cluster
samples according to a predefined grouping factor. It also returns
the table used to map genes to GO terms, the table summarising the
one-way ANOVA results for each gene, and finally the predefined
grouping factor used for ANOVA. Genes annotated to a GO term but
absent from the expression dataset are ignored.

  * get_mart_dataset() returns a connection to the appropriate BioMart
dataset based on the gene name of the first gene in the expression
dataset. The choice of the dataset can be overriden by the user
if a valid BioMart ensembl dataset is specified.
  
  * subset_scores() filters the table of scored GO terms in the
output of GO_anova() and returns a list formatted identically to the 
output of GO_anova() with the resulting filtered table of scores.

  * hist_scores() plots the distribution of average F scores in the
output of GO_anova() or subset_scores().

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
