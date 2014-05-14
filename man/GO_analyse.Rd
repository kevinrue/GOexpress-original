\name{GO_analyse}
\alias{GO_analyse}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Identifies gene ontologies clustering samples according to predefined factor.
}
\description{
Combines local gene expression data with online Gene Ontology (GO) resources to identify ontologies 
enriched for genes best clustering predefined groups of samples based on expression levels.

The package will download Gene Ontologies from the Biomart database for the species
corresponding to the expression dataset prior to the analysis.
A one-way ANOVA is performed to evaluate the ability of each gene to cluster samples according 
to a predefined grouping factor. 
Each GO term is scored and ranked according to the average power of all related genes to
cluster the samples according to the factor of interest.
The ranked list of GO terms is returned, with tools allowing to visualise the statistics
of each gene and ontology.
}
\usage{
GO_analyse(expr_data, phenodata, f, biomart_dataset="", microarray="",
                    adj.P.method = "BH", FDR=0.05)

# To force the use of the human ensembl BioMart, use:
# GO_analyse(expr_data, phenodata, f, biomart_dataset = "hsapiens_gene_ensembl")

# To force use of the bovine affy_bovine microarray annotations use:
# GO_analyse(expr_data, phenodata, f, microarray = "affy_bovine")
}
\arguments{
  \item{expr_data}{
Gene expression matrix. Row names are ensembl gene identifiers or probeset identifiers.
Column names are sample identifiers.
}
  \item{phenodata}{
An AnnotatedDataFrame. Row names are sample idenfifiers in the same order as those in
the expression dataset.
Column names are grouping factors and phenotypic traits usable for the one-way ANOVA.  
}
  \item{f}{
A value in colnames(phenodata) used as the grouping factor for the one-way ANOVA.
}
  \item{biomart_dataset}{
Optional. The BioMart ensembl dataset corresponding to the species studied.
Automatically detected if not specified.
}
  \item{microarray}{
Optional. The BioMart microarray column name in the dataset corresponding to the species studied.
Automatically detected if not specified.
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
A list of four elements containing the results of the analysis.
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
\seealso{
Method \code{\link[biomaRt:getBM]{biomaRt::getBM}}.
}
\examples{
# Valid ensembl datasets are listed in the following variable
data(prefix2dataset)

# Valid microarray= values are listed in the following variable
data(microarray2dataset)


# Load example data subset (100 genes from a RNAseq experiment)
data(log2cpm)

# Load example sample annotation (117 samples)
# Multiple treatments, timepoints, and replicates
data(targets)

# Check that the order of samples are identical in log2cpm and targets
all(colnames(log2cpm) == rownames(targets)) # TRUE

# Run the analysis on factor "Treatment"
# Factor must be one of the column names in pData(targets)
raw_results = GOexpress::GO_analyse(expr_data=log2cpm, phenodata=targets, f="Treatment")
}

\keyword{ GO }
\keyword{ gene }
\keyword{ ontology }
\keyword{expression}
