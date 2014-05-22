\name{expression_plot_symbol}
\alias{expression_plot_symbol}
\title{
Plots the expression profile of a gene by levels of a factor
}
\description{
This function will plot the expression profile of a gene across a valid
X-axis variable from the \code{AnnotatedDataFrame} while representing the mean
and confidence interval of groups of samples defined by levels of a valid
grouping factor from the AnnotatedDataFrame.

In the case of a gene name represented by multiple gene identifiers
in the dataset, a lattice of plots will be produced. Each of the plots in this
lattice can subsequently be plotted separately using its associated index.
}
\usage{
expression_plot_symbol(gene_symbol, expr_data, phenodata, x_var, result,
                                  f=result$factor, index=0, ylab="log2cpm",
                                  col.palette="Accent",
                                  col = RColorBrewer::brewer.pal(n=length(levels(Biobase::pData(phenodata)[,f])), name=col.palette),
                                  level=0.95,
                                  titles=c(), title.size=2)
}
\arguments{
  \item{gene_symbol}{
A gene name present in \code{rownames(raw_result$genes)}.
}
  \item{expr_data}{
A sample-by-gene matrix providing the expression levels of each gene in each sample.
}
  \item{phenodata}{
An \code{AnnotatedDataFrame} from package Biobase providing phenotypic information
about the samples and grouping factors with 2 or more levels.
}
  \item{x_var}{
A value from \code{colnames(pData(phenodata))} to plot on the X-axis. If representing
time on the X-axis, try to store the time-points as numeric values in the
\code{AnnotatedDataFrame} to adequately space the time-points.
}
  \item{result}{
An output of the \code{GO_analyse()} or \code{subset_scores()} function.
}
  \item{f}{
A value from \code{colnames(pData(phenodata))} to group the samples when representing
mean and confidence interval. The factor
specified in the initial \code{GO_analyse()} call is used by default. Unexpected
grouping factors of samples can reveal interesting trends (e.g. "Animal",
"Tissue", "CellType", ...).
}
  \item{index}{
In the case where multiple gene identifiers are associated with the gene
name, \code{index=2} will plot the expression profile of the second gene
identifier alone, for instance.
}
  \item{ylab}{
Change axis labels and legend titles.
}
  \item{col.palette}{
A valid \code{RColorBrewer} palette name to fetch the colormap from.
}
  \item{col}{
A vector of color names or codes. The number of colors provided must match the 
number of levels of the grouping factor. Default to a palette with an adequate
set of colors.
}
  \item{level}{
The confidence interval level to visualise around the mean of each group.
Default is 0.95.
}
  \item{titles}{
Character vector providing as many titles as there are plots to replace
the default titles. Default is a combination of the gene id and the
associated gene.
}
  \item{title.size}{
Changes the font size of the title.
}
}
\value{
  The ggplot object.
}
\references{
  \itemize{
    \item \href{http://cran.r-project.org/web/packages/ggplot2/index.html}{ggplot2} package.
  }
}
\author{
Kevin Rue-Albrecht
}
\section{Warning}{
  Common issues:
    \itemize{
      \item It may not be possible to produce plots where the combination of X-axis
  variable and grouping factor leaves too few replicates to compute
  a confidence interval for each X value. This is a limitation imposed by the ggplot2
  package to produce proper statistics and confidence intervals.
  }
}
\seealso{
  Packages \code{\link[Biobase]{Biobase}} and \code{\link[ggplot2]{ggplot2}}
  , methods
  \code{\link[GOexpress:expression_plot]{GOexpress::expression_plot}} and
  \code{\link[GOexpress:GO_analyse]{GOexpress::GO_analyse}},
  and
  \href{https://github.com/kevinrue/GOexpress/wiki/Typical-analysis-pipeline}{tutorial} on GitHub.
}
\examples{
# load the sample output data
data(raw_result)

# Expression by gene identifier (TNIP3)
GOexpress::expression_plot_symbol(gene_symbol="PIK3AP1", expr_data=log2cpm,
                           phenodata=targets, x_var="Timepoint",
                           result=raw_result)
                           
# Same gene, plotted by animal and grouped by treatment (merging time points)
GOexpress::expression_plot_symbol(gene_symbol="PIK3AP1", expr_data=log2cpm,
                           phenodata=targets, x_var="Animal",
                           result=raw_result,
                           f="Treatment")
                           
# Same gene, plotted by animal and grouped by time-point (merging treatments)
GOexpress::expression_plot_symbol(gene_symbol="PIK3AP1", expr_data=log2cpm,
                           phenodata=targets, x_var="Animal",
                           result=raw_result,
                           f="Time")
}

\keyword{ GOexpress }
\keyword{ gene }
\keyword{ expression }