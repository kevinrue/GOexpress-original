\name{GO_anova}
\alias{GO_anova}
%- Also NEED an '\alias' for EACH other topic documented here.
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
\usage{
# Let the method automatically detect the BioMart dataset to use
GO_anova(expr_data, phenodata, f, biomart_dataset = "", adj.P.method = "BH", FDR = 0.05)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{expr_data}{
Gene expression matrix. Row names are ensembl gene identifiers. Column names are sample identifiers.
}
  \item{phenodata}{
An AnnotatedDataFrame. Row names are sample idenfifiers in the same order as those in the expression dataset.
Column names are grouping factors and phenotypic traits usable for the one-way ANOVA.  
}
  \item{f}{
A value in colnames(phenodata) used as the grouping factor for the one-way ANOVA.
}
  \item{biomart_dataset}{
Optional. The BioMart ensembl dataset corresponding to the species studied. Automatically detected if not specified.
}
  \item{adj.P.method}{
Optional. The multiple testing correction method used to correct the P values obtained in each ANOVA.
}
  \item{FDR}{
Optional. The threshold for significance. Default is 0.05.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A list of two elements containing the results of the analysis.
  \item{ scores }{A table ranking all GO terms related to genes in the expression dataset
  based on the average ability of their related genes to cluster the samples according
  to the predefined grouping factor.}
  \item{ mapping }{The table mapping genes present in the dataset to GO terms.}
  \item{ anova }{A table ranking all genes present in the expression dataset based on 
  their ability to cluster the samples according to the predefined grouping factor}
  \item{ factor }{The predefined grouping factor}
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
%% ~~ Everything here should run (use existing data, or preload some) ~~
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ GO }
\keyword{ gene }
\keyword{ ontology }
%%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
