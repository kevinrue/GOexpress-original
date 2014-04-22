\name{expression_plot}
\alias{expression_plot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plots the expression profile of a gene by levels of a factor
}
\description{
This function will plot the expression profile of a gene across a valid
X-axis variable from the AnnotatedDataFrame while representing the mean
and confidence interval of groups of samples defined by levels of a valid
grouping factor from the AnnotatedDataFrame.
}
\usage{
expression_plot(gene_id, expr_data, phenodata, x_var, result, f=result$factor, 
                           ylab = "log2(cpm)", col.palette="Accent",
                           col = brewer.pal(n=length(levels(pData(phenodata)[,f])), name=col.palette),
                           level=0.95, title=NULL, title.size=2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gene_id}{
An gene identifier present in rownames(expr_data).
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
An output of the \code{GO_anova()} or \code{subset_scores()} function.
}
  \item{f}{
A value from colnames(pData(phenodata)) to group the samples when representing
mean and confidence interval. If the \code{result} argument is used, the factor
specified in the initial \code{GO_anova()} call is used by default. Unexpected
grouping factors of samples can reveal interesting trends (e.g. "Animal",
"Tissue", "CellType").
}
  \item{ylab}{
Change axis labels and legend titles.
}
  \item{col.palette}{
A valid RColorBrewer palette name to fetch the colormap from.
}
  \item{col}{
A vector of color names or codes. The number of colors provided must match the 
number of levels of the grouping factor. Default to a color-blind friendly
palette of three colors.
}
  \item{level}{
The confidence interval level to visualise around the mean of each group.
Default is 0.95.
}
  \item{title}{
Changes the plot title. Default is a combination of the gene id and the
associated gene.
}
  \item{title.size}{
Changes the font size of the title.
}
}
\details{
If the \code{result} argument is not specified and no custom title is specified either,
the function will connect to the BioMart dataset for the species corresponding
to the gene identifier to generate a title combining the gene identifier and 
the associated gene name.
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

