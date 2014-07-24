\name{cluster_GO}
\alias{cluster_GO}
\title{
Generates a hierarchical clustering / dendrogram of the samples
}
\description{
Clusters the samples using only the expression levels of genes associated with a
given go_id.
}
\usage{
cluster_GO(go_id, result, expr_data, phenodata, f=result$factor, 
           method_dist="euclidean", method_hclust="average", cex=0.8,
           main=paste(go_id, result$GO[result$GO$go_id == go_id, "name_1006"]),
           xlab="Distance",
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
Gene-by-sample expression matrix. Row names are ensembl gene identifiers or probeset identifiers.
Column names are sample identifiers.
}
  \item{phenodata}{
An \code{AnnotatedDataFrame}. Row names are sample idenfifiers in the same order as those in the expression dataset.
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
\value{
  Returns the output of the \code{plot()} function.
}
\author{
Kevin Rue-Albrecht
}

\seealso{
  Method \code{\link[GOexpress:GO_analyse]{GO_analyse}}
  or
  \href{https://github.com/kevinrue/GOexpress/wiki/Typical-analysis-pipeline}{tutorial}
  on GitHub.
}
\examples{
# load the sample output data
data(raw_result)

# Hierarchical clustering of samples based on the same GO term
cluster_GO(go_id="GO:0034142", result=raw_result, expr_data=log2cpm,
                      phenodata=targets, cex=0.7)

# Re-label sample by another factor
cluster_GO(go_id="GO:0034142", result=raw_result, expr_data=log2cpm,
                      phenodata=targets, cex=0.7,
                      f="Group")

}

\keyword{ GOexpress }
\keyword{ gene }
\keyword{ expression }
\keyword{ clustering }
