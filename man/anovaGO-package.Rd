\name{GOexpress-package}
\alias{GOexpress-package}
\alias{GOexpress}
\docType{package}
\title{
Visualise microarray and RNAseq data with gene ontology annotations.
}
\description{
Combines local gene expression data with online Gene Ontology (GO)
resources to visualise and rank ontologies enriched for genes best clustering
predefined groups of samples based on gene expression levels. The
package will import Gene Ontologies from the Biomart database for
the species corresponding to the expression dataset prior to
performing the analysis. A randomForest analysis is performed to evaluate
the ability of each gene to cluster samples according to a predefined
grouping factor (one-way ANOVA available as an alteranative). Each GO
term is scored and ranked according to the average power of all associated
genes to cluster the samples according to the factor. The ranked list
of GO terms is returned, with tools allowing to visualise the statistics
on a gene- and ontology-basis.
}
\details{
 \tabular{ll}{
 Package: \tab GOexpress\cr
 Type: \tab Package\cr
 Version: \tab 0.99.0\cr
 Date: \tab 2014-05-22\cr
 License: \tab  GPL (>= 3)\cr
 }
 This package requires only two input variables
 \enumerate{
  \item A sample-by-gene matrix providing gene expression values for each gene in each sample
  \item An AnnotatedDataFrame from the Biobase package providing phenotypic information and
 grouping factors with two or more levels.
 }
 Following analysis, visualisation methods include:
 \itemize{
  \item Histogram and quantiles representations of the scores of GO terms
  \item Filtering of results on various criteria (e.g. number of genes in GO term)
  \item Re-ordering of GO terms and gene result tables based on score or rank metric
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
\examples{
  # Sample input data available with package:
  data(log2cpm)
  data(targets)
  
  # Sample output data available with package:
  data(raw_result)
  
  # Supported species and microarrays:
  data(microarray2dataset)
  data(prefix2dataset)
}
\references{
\itemize{
  \item \href{http://cran.r-project.org/web/packages/gplots/index.html}{gplots} package.
  \item \href{http://cran.r-project.org/web/packages/ggplot2/index.html}{ggplot2} package.
  \item \href{http://cran.r-project.org/web/packages/randomForest/index.html}{randomForest} package.
  \item \href{http://cran.r-project.org/web/packages/stringr/index.html}{stringr} package.
  \item \href{http://cran.r-project.org/web/packages/VennDiagram/index.html}{VennDiagram} package.
  }
}

\keyword{ package }
\keyword{ GOexpress }
\keyword{ randomForest }
\keyword{ anova }
\keyword{ gene }
\keyword{ expression }
\keyword{ clustering }
\keyword{ ontology }

\seealso{
  Main method for an example usage:
  \code{\link[GOexpress:GO_analyse]{GOexpress::GO_analyse}}
  or
  \href{https://github.com/kevinrue/GOexpress/wiki/Typical-analysis-pipeline}{tutorial} on GitHub.

  Packages
  \code{\link[Biobase]{Biobase}},
  \code{\link[ggplot2]{ggplot2}},
  \code{\link[randomForest]{randomForest}},
  \code{\link[RColorBrewer]{RColorBrewer}},
  \code{\link[VennDiagram]{VennDiagram}}.
  
  Methods
  \code{\link[biomaRt:getBM]{biomaRt::getBM}},
  \code{\link[gplots:heatmap.2]{gplots::heatmap.2}},
  \code{\link[stringr:str_extract]{stringr::str_extract}},.
}
