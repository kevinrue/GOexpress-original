\name{heatmap_GO}
\alias{heatmap_GO}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generates a heatmap and hierarchical clustering of the samples and the genes
}
\description{
Clusters the samples and the genes associated with a GO term using the expression 
levels of genes related to a given ontology. Represents expression levels of those genes
in a heatmap.
}
\usage{
heatmap_GO(go_id, result, expr_data, phenodata, gene_names=TRUE,
           f=result$factor, scale="none", cexCol=1.2, cexRow=0.5, 
           trace="none", expr.col=rev(redgreen(75)), row.col.palette="Accent",
           row.col=brewer.pal(n=length(unique(pData(phenodata)[,f])), name=row.col.palette),
           main=paste(go_id, result$scores[result$scores$go_id == go_id,]$name_1006),
           ...)
}
\arguments{
  \item{go_id}{
A Gene Ontology (GO) identifier.
}
  \item{result}{
The output of \code{GO_analyse()} or a subset of it obtained from \code{subset_scores()}.
}
  \item{expr_data}{
Gene expression matrix. Row names are ensembl gene identifiers or probeset identifiers.
Column names are sample identifiers.
}
  \item{phenodata}{
An AnnotatedDataFrame. Row names are sample idenfifiers in the same order as those in
the expression dataset.
Column names are grouping factors and phenotypic traits usable for the one-way ANOVA.  
}
  \item{gene_names}{
A boolean value. Default is TRUE, to label genes by their associated gene name.
If FALSE, labels the genes by their gene identifier in the expression dataset.
}
  \item{f}{
The predefined grouping factor.
}
  \item{scale}{
Character indicating if the values should be centered and scaled in either the row
direction or the column direction, or none. Default is "none". See \code{heatmap.2()}.
}
  \item{cexCol, cexRow}{
Positive numbers, used as cex.axis in for the row or column axis labeling. Defaults are 
.2 and 1, respectively. See \code{heatmap.2()}.
}
  \item{trace}{
Character string indicating whether a solid "trace" line should be drawn across 'row's
or down 'column's, both' or 'none'. The distance of the line from the center of each
color-cell is proportional to the size of the measurement. Defaults to 'none'
}
  \item{expr.col}{
Character vector indicating the colors to represent the different levels of gene
expression. Defaults to a colormap of 75 shades ranging from green (low) to red
(high) and centered around black.
}
  \item{row.col.palette}{
A valid RColorBrewer palette name to fetch the colormap from, to color-code the groups
of samples.
}
  \item{row.col}{
A vector of color names or codes. The number of colors provided must match the 
number of levels of the grouping factor. Default to an adequate palette of colors.
}
  \item{main}{
Main title of the figure. Default is paste(go_id, go_name).
}
  \item{\dots}{
Additional arguments passed on to \code{heatmap.2()}.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  Returns the output of the \code{gplots::heatmap.2()} function.
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
\seealso{
Method \code{\link[gplots:heatmap.2]{gplots::heatmap.2}}.
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
