\name{table_genes}
\alias{table_genes}
\title{
Returns a table listing the genes associated with a given Gene Ontology
}
\description{
Given a Gene Ontology (GO) identifier represented in the dataset and the output
variable of a \code{GO_analyse()} function, \code{table_genes()} returns a
table listing the genes associated with that go_id, their associated gene name,
and description.
}
\usage{
table_genes(go_id, result, data.only=FALSE)
}
\arguments{
  \item{go_id}{
A Gene Ontology (GO) identifier.
}
  \item{result}{
The output of the \code{GO_analyse()} function.
}
  \item{data.only}{
Whether to return only the gene identifiers present in the given dataset
or alternatively returns all gene identifiers associated with the GO term
in the BioMart database.
}
}
\value{
A data frame listing the statistics and BioMart annotations for the genes 
present in the expression dataset and associated with the GO term.
}
\author{
Kevin Rue-Albrecht
}
\seealso{
  Method \code{\link[GOexpress:GO_analyse]{GOexpress::GO_analyse}}.
}
\examples{
# load the sample output data
data(raw_result)

# Table of result for genes associated with the GO term
# "toll-like receptor 4 signaling pathway"
GOexpress::table_genes(result=raw_result, go_id="GO:0034142")
}

\keyword{ GOexpress }
\keyword{ gene }
\keyword{ ontology }
