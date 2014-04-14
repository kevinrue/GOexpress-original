\name{cluster_GO}
\alias{cluster_GO}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generates a hierarchical clustering dendrogram of the samples
}
\description{
Clusters the samples using the expression levels of genes related to a
given go_id.
}
\usage{
cluster_GO(go_id, result, expr_data, phenodata, f=result$factor, 
                      method_dist="euclidean", method_hclust="average", cex=0.8,
                      main=paste(go_id, result$scores[result$scores$go_id == go_id,]$name_1006),
                      xlab="Distance", ...)

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
  \item{f}{
The predefined grouping factor.
}
  \item{method_dist}{
The method used to calculate distance between samples. See the 
\code{dist()} method from package \code{stats}.
}
  \item{method_hclust}{
The method used to cluster samples. See the 
\code{hclust()} method from package \code{stats}.
}
  \item{cex}{
A numeric value defining the character expansion of text in the plot.
}
  \item{main}{
A character string for the main title of the plot.
}
  \item{xlab}{
A label for the x axis, defaults to "Distance".
}
  \item{\dots}{
Additional parameters passed on to \code{dist()}, \code{hclust()}
and \code{plot()}.
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
\code{\link[anovaGO:GO_anova]{anovaGO::GO_anova}}
}
\examples{

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ GO }
\keyword{ gene }
\keyword{ ontology }
%%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
