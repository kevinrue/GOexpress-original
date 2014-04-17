\name{expression_plot_symbol}
\alias{expression_plot_symbol}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plots the expression profile of a gene by levels of a factor
}
\description{
This function will plot the expression profile of a gene across a valid
X-axis variable from the AnnotatedDataFrame while representing the mean
and confidence interval of groups of samples defined by levels of a valid
grouping factor from the AnnotatedDataFrame.

In the case of a gene name represented by multiple ensembl gene identifiers
in the dataset, a lattice of plots will be produced. Each of the plots in this
lattice can be plotted separately using its associated index.
}
\usage{
expression_plot_symbol(gene_symbol, expr_data, phenodata, x_var,
                       result = NULL, f = result$factor, index = 0,
                       biomart_dataset = "", ylab="log2cpm",
                       cbPalette = c("#56B4E9", "#D55E00", "#F0E442"),
                       level = 0.95, titles = c(), title.size = 2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gene_symbol}{
A gene name present in rownames(result$anova) if results is specified.
Otherwise, the BioMart database will be queried for ensembl gene
identifiers associated with the gene name.
}
  \item{expr_data}{
A sample-by-gene matrix providing expression levels of each gene in each sample.
}
  \item{phenodata}{
An AnnotatedDataFrame from package Biobase providing phenotypic information
about the samples and grouping factors with 2 or more levels.
}
  \item{x_var}{
A value from colnames(pData(phenodata)) to plot on the X-axis. If representing
time on the X-axis, make sure the time-points are stored as numeric values in the
AnnotatedDataFrame for best results.
}
  \item{result}{
An output of the \code{GO_anova()} or \code{subset_scores()} function. If not
specified, the BioMart database will be queried to look for ensembl gene 
identifiers associated with the gene name, and the \code{f=} argument is
then mandatory.
}
  \item{f}{
Optional. A value from colnames(pData(phenodata)) to group the samples when representing
mean and confidence interval. If the \code{result} argument is used, the factor
specified in the initial \code{GO_anova()} call is used by default. Unexpected
grouping factors of samples can reveal interesting trends (e.g. "Animal",
"Tissue", "CellType").
}
  \item{index}{
In the case where multiple ensembl gene identifiers are associated with the gene
name, \code{index=2} will plot the expression profile of the second ensembl
identifier alone, for instance.
}
  \item{biomart_dataset}{
Optional. The BioMart ensembl dataset corresponding to the species studied.
Automatically detected if not specified.
}
  \item{ylab}{
Change axis labels and legend titles.
}
  \item{cbPalette}{
Optional. A vector of color names or codes. The number of colors provided must match the 
number of levels of the grouping factor. Default to a color-blind friendly
palette of three colors.
}
  \item{level}{
The confidence interval level to visualise around the mean of each group.
Default is 0.95.
}
  \item{titles}{
Character vector providing as many titles as there are plots to replace
the default titles. Default is a combination of the ensembl gene id and the
associated gene.
}
  \item{title.size}{
Changes the font size of the title.
}
}
\details{
If provided, the \code{result=} argument must contain the output of either the 
\code{GO_anova()} or the \code{subset_scores()} functions. Those objects contain
useful information such as the mapping between gene name and corresponding  
ensembl gene identifiers, but also the factor used for the \code{GO_anova()} function.
This avoids the need to fetch gene annotations from BioMart and also provides a
default grouping factor for the expression profile.

However, if the \code{result=} argument is not provided, gene annotations will
be queried from BioMart, and a grouping factor must be provided by the user.

}
\value{
  The ggplot object.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Kevin Rue-Albrecht
}
\note{
%%  ~~further notes~~
}
\section{Warning}{
  Common issues:
    \itemize{
      \item It may not be possible to produce plots where the combination of X-axis
  variable and grouping factor leaves less than 10 replicates to compute
  a confidence interval for each X value. This is a limitation imposed by the ggplot2
  package to produce proper statistics and confidence intervals.
  }
}
\seealso{
Packages \code{\link[Biobase]{Biobase}} and \code{\link[ggplot2]{ggplot2}}.
}
\examples{
%%##---- Should be DIRECTLY executable !! ----
%%##-- ==>  Define data, use random,
%%##--  or do  help(data=index)  for the standard data sets.
}

\keyword{ GO }
\keyword{ gene }
\keyword{ ontology }
\keyword{expression}
