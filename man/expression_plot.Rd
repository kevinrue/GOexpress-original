\name{expression_plot}
\alias{expression_plot}
\title{
Plots the expression profile of a gene by levels of a factor
}
\description{
This function will plot the expression profile of a gene across a valid
X-axis variable from the \code{AnnotatedDataFrame} while representing the mean
and confidence interval of groups of samples defined by levels of a valid
grouping factor from the AnnotatedDataFrame.
}
\usage{
expression_plot(
    gene_id, result, eSet, x_var, 
    f=result$factor, ylab = "log2(cpm)", col.palette="Accent",
    col = brewer.pal(n=length(levels(pData(eSet)[,f])), name=col.palette),
    level=0.95, title=NULL, title.size=2, axis.title.size=20,
    axis.text.size=15, legend.text.size=15, legend.title.size=20,
    legend.key.size=30)
}
\arguments{
  \item{gene_id}{
An gene or probeset identifier present in \code{rownames(expr_data)}.
}
  \item{result}{
An output of the \code{GO_analyse()} or \code{subset_scores()} function.
}
  \item{eSet}{
\code{ExpressionSet} of the \code{Biobase} package including a gene-by-sample
expression matrix in the \code{AssayData} slot, and phenotypic information
data-frame in the \code{phenodate} slot. In the expression matrix, row names
are Ensembl gene identifiers or probeset identifiers, and column names are
sample identifiers. In the phentypic data-frame, row names are sample
idenfifiers, column names are grouping factors and phenotypic traits usable
for the one-way ANOVA. 
}
  \item{x_var}{
A value from \code{colnames(pData(phenodata))} to plot on the X-axis. If
representing
time on the X-axis, try to store the time-points as numeric values in the
\code{AnnotatedDataFrame} to adequately space the time-points.
}
  \item{f}{
A value from \code{colnames(pData(phenodata))} to group the samples when
representing
mean and confidence interval. The factor
specified in the initial \code{GO_analyse()} call is used by default.
Unexpected grouping factors of samples can reveal interesting trends (e.g.
"Animal", "Tissue", "CellType", ...).
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
  \item{title}{
Changes the plot title. Default is a combination of the gene id and the
associated gene.
}
  \item{title.size}{
Changes the font size of the title.
}
  \item{axis.title.size}{
Changes the font size of the axes title.
}
  \item{axis.text.size}{
Changes the font size of the axes text labels.
}
  \item{legend.text.size}{
Changes the font size of the legend text labels.
}
  \item{legend.title.size}{
Changes the font size of the legend title.
}
  \item{legend.key.size}{
Changes the size of the legend keys.
}
}
\details{
If the \code{result} argument is not specified and no custom title is
specified either, the function will combine the feature identifier and the
associated gene name.
}
\value{
  The ggplot object.
}
\references{
  \itemize{
    \item \href{http://cran.r-project.org/web/packages/ggplot2/index.html}{
    ggplot2} package.
  }
}
\author{
Kevin Rue-Albrecht
}
\section{Warning}{
  Common issues:
    \itemize{
      \item It may not be possible to produce plots where the combination of
      X-axis variable and grouping factor leaves too few replicates to compute
  a confidence interval for each X value. This is a limitation imposed by the
  ggplot2 package to produce proper statistics and confidence intervals.
  }
}
\seealso{
  Packages \code{\link[Biobase]{Biobase}} and \code{\link[ggplot2]{ggplot2}}
  , methods
  \code{\link[GOexpress:expression_plot_symbol]{
  expression_plot_symbol}}
  and
  \code{\link[GOexpress:expression_plot_symbol]{
  expression_plot_symbol}}.
}
\examples{
# load the sample output data
data(raw_results)

# Expression by gene identifier (TNIP3)
expression_plot(gene_id="ENSBTAG00000047107", eSet=AlvMac, x_var="Timepoint",
                           result=raw_results)
                           
# Same gene, plotted by animal and grouped by treatment (merging time points)
expression_plot(gene_id="ENSBTAG00000047107", eSet=AlvMac, x_var="Animal",
                           result=raw_results,
                           f="Treatment")
                           
# Same gene, plotted by animal and grouped by time-point (merging treatments)
expression_plot(gene_id="ENSBTAG00000047107", eSet=AlvMac, x_var="Animal",
                           result=raw_results,
                           f="Time")
}

\keyword{ GOexpress }
\keyword{ gene }
\keyword{ expression }
