\name{anovaGO-package}
\alias{anovaGO-package}
\alias{anovaGO}
\docType{package}
\title{
Identifies gene ontologies clustering samples according to predefined factor.
}
\description{
Combines local gene expression data with online Gene Ontology (GO) resources to identify ontologies 
enriched for genes best clustering predefined groups of samples based on expression levels.
The package will download Gene Ontologies from the Biomart database for the species
corresponding to the expression dataset priori to performing the analysis.
A one-way ANOVA is performed to evaluate the ability of each gene to cluster samples according 
to a predefined grouping factor. 
Each GO term is scored and ranked according to the average power of all related genes to
cluster the samples according to the factor.
The ranked list of GO terms is returned, with tools allowing to visualise the statistics
on a gene- and ontology-bases.
}
\details{
 \tabular{ll}{
 Package: \tab anovaGO\cr
 Type: \tab Package\cr
 Version: \tab 0.3.1\cr
 Date: \tab 2014-04-16\cr
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
  \item Table of statistics for genes annotated to a GO term
  \item Hierarchical clustering of samples based on the expression level of genes
  annotated to a GO term
  \item Heatmap of samples and genes based on the expression level of genes
  annotated to a GO term
  \item Expression profile of a gene against one factor (e.g. Time) while grouping 
  samples on another factor (e.g. Treatment)
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
  See also: package \code{\link[ggplot2]{ggplot2}} and method
  \code{\link[biomaRt:getBM]{biomaRt::getBM}}, .
}
\examples{
%%##---- Should be DIRECTLY executable !! ----
%%##-- ==>  Define data, use random,
%%##--  or do  help(data=index)  for the standard data sets.
}
