\name{GOexpress-package}
\alias{GOexpress-package}
\alias{GOexpress}
\docType{package}
\title{
Visualise microarray and RNAseq data with gene ontology annotations.
}
\description{
Combines local gene expression data with online Gene Ontology (GO) resources to identify ontologies 
enriched for genes best clustering predefined groups of samples based on gene expression levels.
The package will import Gene Ontologies from the Biomart database for the species
corresponding to the expression dataset prior to performing the analysis.
A one-way ANOVA is performed to evaluate the ability of each gene to cluster samples according 
to a predefined grouping factor (alternate metrics currently being tested). 
Each GO term is scored and ranked according to the average power of all related genes to
cluster the samples according to the factor.
The ranked list of GO terms is returned, with tools allowing to visualise the statistics
on a gene- and ontology-basis.
}
\details{
 \tabular{ll}{
 Package: \tab GOexpress\cr
 Type: \tab Package\cr
 Version: \tab 0.5.5\cr
 Date: \tab 2014-05-15\cr
 License: \tab GPL-2\cr
 }
 This package requires only two input variables
 \enumerate{
  \item A sample-by-gene matrix providing gene expression values for each gene in each sample
  \item An AnnotatedDataFrame from the Biobase package providing phenotypic information and
 grouping factors with two or more levels.
 }
 Following analysis, visualisation tools include:
 \itemize{
  \item Histogram and quantiles representations of the scores of GO terms
  \item Filtering of results on various criteria (e.g. number of genes in GO term)
  \item Table of statistics for genes annotated to a given GO term
  \item Hierarchical clustering of samples based on the expression level of genes
  annotated to a given GO term
  \item Heatmap of samples and genes based on the expression level of genes
  annotated to a given GO term
  \item Expression profile of a gene against one given factor (e.g. Time) while grouping 
  samples on another given factor (e.g. Treatment)
  \item Univariate analysis of the expression level of a gene in the different groups
  of each experimental factor.
  \item Venn diagram of the genes shared between a list of GO terms.
 }
}
\author{
Kevin Rue-Albrecht

Maintainer: Kevin Rue-Albrecht <kevin.rue@ucdconnect.ie>
}
\references{
%% ~put references to the literature/web site here ~
}
\keyword{ package }
\seealso{
  See also: packages
  \code{\link[Biobase]{Biobase}},
  \code{\link[ggplot2]{ggplot2}},
  
  \code{\link[RColorBrewer]{RColorBrewer}},
  \code{\link[VennDiagram]{VennDiagram}},
  and methods
  \code{\link[biomaRt:getBM]{biomaRt::getBM}},
  \code{\link[gplots:heatmap.2]{gplots::heatmap.2}}.
}
\examples{
%%##---- Should be DIRECTLY executable !! ----
%%##-- ==>  Define data, use random,
%%##--  or do  help(data=index)  for the standard data sets.
}
