\name{list_genes}
\alias{list_genes}
\title{
Returns the genes associated with a Gene Ontology
}
\description{
Given a Gene Ontology (GO) identifier represented in the dataset,
returns a character vector listing the gene identifiers annotated to it.
}
\usage{
list_genes(go_id, result, data.only=TRUE)
}
\arguments{
  \item{go_id}{
A Gene Ontology (GO) identifier represented in the dataset.
}
\item{result}{
The output of \code{GO_analyse()} or a subset of it obtained from
\code{subset_scores()}.
}
  \item{data.only}{
Whether to return only the gene identifiers present in the given dataset
or alternatively returns all gene identifiers associated with the GO term
in the BioMart database. Default is TRUE.
}
}
\value{
A character vector listing the gene identifiers of the genes associated
with the GO term.
}
\author{
Kevin Rue-Albrecht
}
\seealso{
  Method \code{\link[GOexpress:GO_analyse]{GO_analyse}} and
  \href{https://github.com/kevinrue/GOexpress/wiki/Typical-analysis-pipeline}{Tutorial}
  on GitHub.
}
\examples{
# load the sample output data
data(raw_result)

# List of genes associated with the GO term "toll-like receptor 4 signaling pathway"
list_genes(result=raw_result, go_id="GO:0034142")
}

\keyword{ GOexpress }
\keyword{ gene }