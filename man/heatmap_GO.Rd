\name{heatmap_GO}
\alias{heatmap_GO}
\title{
Generates a heatmap and hierarchical clustering of the samples and the genes
}
\description{
Clusters the samples and the genes associated with a GO term using the
expression levels of genes related to a given ontology. Represents expression
levels of those genes in a heatmap.
}
\usage{
heatmap_GO(
    go_id, result, expr_data, phenodata, gene_names=TRUE,
    f=result$factor, scale="none", cexCol=1.2, cexRow=0.5, trace="none",
    expr.col=bluered(75), row.col.palette="Accent",
    row.col=brewer.pal(n=length(unique(pData(phenodata)[,f])),
                       name=row.col.palette),
    main=paste(go_id, result$GO[result$GO$go_id == go_id, "name_1006"]), ...)
}
\arguments{
  \item{go_id}{
A Gene Ontology (GO) identifier.
}
  \item{result}{
The output of \code{GO_analyse()} or a subset of it obtained from
\code{subset_scores()}.
}
  \item{expr_data}{
Gene expression matrix. Row names are Ensembl gene identifiers or probeset
identifiers. Column names are sample identifiers.
}
  \item{phenodata}{
An AnnotatedDataFrame. Row names are sample idenfifiers in the same order as
those in the expression dataset. Column names are grouping factors and
phenotypic traits used for analysis and plotting.  
}
  \item{gene_names}{
A boolean value. Default is TRUE, to label genes by their associated gene name.
If FALSE, labels the genes by their feature identifier in the expression
dataset (i.e. ensembl gene identifier or microarray probeset).
}
  \item{f}{
The predefined grouping factor.
}
  \item{scale}{
Character indicating if the values should be centered and scaled in either the
row direction or the column direction, or none. Default is "none".
See \code{heatmap.2()}.
}
  \item{cexCol, cexRow}{
Positive numbers, used as cex.axis in for the row or column axis labeling.
Defaults are 1.2 and 1, respectively. See \code{heatmap.2()}.
}
  \item{trace}{
Character string indicating whether a solid "trace" line should be drawn across
each 'row' or down each 'column', both' or 'none'. The distance of the line
from the center of each color-cell is proportional to the size of the
measurement. Defaults to 'none'.
}
  \item{expr.col}{
Character vector indicating the colors to represent the different levels of
gene expression. Defaults to a colormap of 75 shades ranging from blue (low) to
red (high) and centered around white. If using differential expression data,
you should probably use \code{greenred(75)} instead.
}
  \item{row.col.palette}{
A valid \code{RColorBrewer} palette name to fetch the colormap from, to
color-code the groups of samples.
}
  \item{row.col}{
A vector of color names or codes. The number of colors provided must match the 
number of levels of the grouping factor. Default to an palette of up to 9
colors marking the different levels of the predefined grouping factor on the
left side of the heatmap.
}
  \item{main}{
Main title of the figure. Default is paste(go_id, go_name).
}
  \item{\dots}{
Additional arguments passed on to \code{heatmap.2()}.
}
}
\value{
  Returns the output of the \code{heatmap.2()} function.
}
\references{
\itemize{
  \item \href{http://cran.r-project.org/web/packages/gplots/index.html}{gplots}
  package.
  }
}
\author{
Kevin Rue-Albrecht
}
\seealso{
  Method \code{\link[gplots:heatmap.2]{heatmap.2}},
  \code{\link[GOexpress:GO_analyse]{GO_analyse}},
  \code{\link[RColorBrewer:brewer.pal.info]{brewer.pal.info}}
  and
  \href{https://github.com/kevinrue/GOexpress/wiki/Typical-analysis-pipeline}{tutorial}
  on GitHub.
}
\examples{
# load the sample output data
data(raw_result)

# Heatmap the top-ranked GO term (toll-like receptor 4 signaling pathway) as
# example
heatmap_GO("GO:0034142", raw_result, log2cpm, targets)

# Same with larger sample labels on the right hand side.
heatmap_GO("GO:0034142", raw_result, log2cpm, targets, cexRow=1)

# Change the color-coding to green-black-red gradient (more appropriate for
# differential expression values)
library(gplots)
heatmap_GO("GO:0034142", raw_result, log2cpm, targets, expr.col=greenred(75))
}

\keyword{ GOexpress }
\keyword{ gene }
\keyword{ expression }
\keyword{ clustering }
\keyword{ ontology }