\name{heatmap_GO}
\alias{heatmap_GO}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generates a heatmap and hierarchical clustering of the samples and the genes
}
\description{
Clusters the samples and the genes associated with a GO term using the expression 
levels of genes related to a given go_id. Represents expression levels of those genes
in a heatmap.
}
\usage{
heatmap_GO(go_id, result, expr_data, phenodata, gene_names = TRUE, f = result$factor, scale = "none", cexCol = 1.2, cexRow = 1, main = paste(go_id, result$scores[result$scores$go_id == go_id, ]$name_1006), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{go_id}{
A Gene Ontology (GO) identifier.
}
  \item{result}{
The output of \code{GO_anova()} or a subset of it obtained from \code{subset_scores()}.
}
  \item{expr_data}{
Gene expression matrix. Row names are ensembl gene identifiers. Column names are sample identifiers.
}
  \item{phenodata}{
An AnnotatedDataFrame. Row names are sample idenfifiers in the same order as those in the expression dataset.
Column names are grouping factors and phenotypic traits usable for the one-way ANOVA.  
}
  \item{gene_names}{
A boolean value. Default is TRUE, to label genes by their associated gene name.
If FALSE, labels the genes by their ensembl_gene_id.
}
  \item{f}{
The predefined grouping factor.
}
  \item{scale}{
character indicating if the values should be centered and scaled in either the row direction or the column direction, or none. Default is "none". See \code{heatmap()}.
}
  \item{cexCol, cexRow}{
positive numbers, used as cex.axis in for the row or column axis labeling. Defaults are 
1.2 and 1, respectively. See \code{heatmap()}.
}
  \item{main}{
Main title of the figure. Default is paste(go_id, go_name).
}
  \item{\dots}{
Additional arguments passed on to \code{heatmap()}.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
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

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ GO }
\keyword{ gene }
\keyword{ ontology }
%%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
