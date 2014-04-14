\name{anovaGO-package}
\alias{anovaGO-package}
\alias{anovaGO}
\docType{package}
\title{
Identifies gene ontologies clustering samples according to predefined factor.
}
\description{
Combines gene expression data and Gene Ontology (GO) annotations to identify ontologies 
enriched for genes best clustering predefined groups of samples base on expression levels.
Genes annotated to an ontology (Biological Process) are obtained from the BioMart database.
A one-way ANOVA is performed to evaluate the ability of each gene to cluster samples according 
to a predefined grouping factor.
Each GO term is scored and ranked according to the average power of all related genes to
cluster the samples according to the factor.
The ranked list of GO terms is returned, with tools allowing to visualise the statistics
on a gene- and ontology-bases.
}
\details{
\tabular{ll}{
Package: \tab anovaGO\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2014-04-13\cr
License: \tab GPL-v2\cr
}
An overview of how to use the package, including the most important functions
}
\author{
Kevin Rue-Albrecht

Maintainer: Kevin Rue-Albrecht <kevin.rue@ucdconnect.ie>
}
\references{

}

\keyword{ package }
\seealso{
See also:\code{\link[biomaRt:getBM]{biomaRt:getBM}}
}
\examples{

}
