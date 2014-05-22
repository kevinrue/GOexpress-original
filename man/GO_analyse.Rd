\name{GO_analyse}
\alias{GO_analyse}
\title{
Identifies gene ontologies clustering samples according to predefined factor.
}
\description{
Combines local gene expression data with online Gene Ontology (GO) resources to identify ontologies 
enriched for genes best clustering predefined groups of samples based on expression levels.

The package will download Gene Ontologies from the Biomart database for the species
corresponding to the expression dataset prior to the analysis.
A random forest analysis is performed to evaluate the ability of each gene to
cluster samples according  to a predefined grouping factor. 
Each GO term is scored and ranked according to the average power of all associated genes to
cluster the samples according to the factor of interest. The ranked list of GO terms
is returned, with tools allowing to visualise the statistics
of each gene and ontology.
}
\usage{
GO_analyse(expr_data, phenodata, f, biomart_dataset="", microarray="",
           method="randomForest", rank.by="rank",
           do.trace=100, ntree=1000, mtry=ceiling(2*sqrt(nrow(expr_data))), 
           ...)

# To force the use of the human ensembl BioMart, use:
# GO_analyse(expr_data, phenodata, f, biomart_dataset = "hsapiens_gene_ensembl")

# To force use of the bovine affy_bovine microarray annotations use:
# GO_analyse(expr_data, phenodata, f, microarray = "affy_bovine")
}
\arguments{
  \item{expr_data}{
Gene expression matrix. Row names are ensembl gene identifiers or probeset identifiers.
Column names are sample identifiers. Typically the output of \code{edgeR::cpm()} or
\code{Biobase::exprs()}
}
  \item{phenodata}{
An \code{AnnotatedDataFrame}. Row names are sample idenfifiers in the same order as those in
the expression dataset.
Column names are grouping factors and phenotypic traits used for analysis and plotting.  
}
  \item{f}{
A value in colnames(phenodata) used as the grouping factor for the analysis.
}
  \item{biomart_dataset}{
The BioMart ensembl dataset corresponding to the species studied.
Automatically detected if not specified. Use data(prefix2dataset) for valid values.
}
  \item{microarray}{
The BioMart microarray column name in the dataset corresponding to the species studied.
Automatically detected if not specified. Use data(microarray2dataset) for valid values.
}
  \item{method}{
Either "randomForest" or "rf" to use the random forest framework, or alternatively
either of "anova" or "a" to use the one-way ANOVA model. Default is "randomForest".
}
  \item{rank.by}{
Either of "rank" or "score" to chose the metric used to order the gene and GO term
result tables.
}
  \item{do.trace}{
Only used if method="randomForest". If set to TRUE, give a more verbose output as
randomForest is run. If set to some integer, then running output is printed for every
do.trace trees.
}
  \item{ntree}{
Only used if method="randomForest". Number of trees to grow. This should not be set to
too small a number, to ensure that every input row gets predicted at least a few times
}
  \item{mtry}{
Only used if method="randomForest". Number of variables randomly sampled as candidates
at each split. Default value is 2*sqrt(gene_count) which is approximately 220 genes
for a dataset of 12,000 genes.
}
  \item{\dots}{
Additional arguments passed on to the randomForest() method.
}
}
\value{
A list containing the results of the analysis. Some elements only appear in the output
of a randomForest analysis. Core elements:
  \item{ GO }{A table ranking all GO terms related to genes in the expression dataset
  based on the average ability of their related genes to cluster the samples according
  to the predefined grouping factor.}
  \item{ mapping }{The table mapping genes present in the dataset to GO terms.}
  \item{ genes }{A table ranking all genes present in the expression dataset based on 
  their ability to cluster the samples according to the predefined grouping factor.}
  \item{ factor }{The predefined grouping factor.}
  \item{ method }{The statistical framework used.}
}
\details{
  The current scoring functions strongly favor GO terms associated with fewer genes at 
  the top of the ranking. This bias may actually be seen as a valuable feature which
  enables the user to browse through GO terms of increasing "granularity", i.e.
  associated with increasingly large sets of genes, although consequently being
  increasingly vague and blurry (e.g. "protein binding" molecular function associated
  with over 6,000 genes).
  
  It is suggested to use the \code{GOexpress::subset_scores()} function to filter out
  GO terms with fewer than 5+ genes associated with it. Indeed, those GO terms are
  more sensitive to outlier genes as they were scored on the average of a handful of
  genes.
}
\references{
\itemize{
  \item \href{http://cran.r-project.org/web/packages/randomForest/index.html}{randomForest} package.
  }
}
\author{
Kevin Rue-Albrecht
}
\section{Warning}{
  Make sure that the factor \code{f} is an actual factor is the R language meaning.
  This is important for the underlying statistical framework to identify the groups
  of samples defined by their level of this factor.
  
  If the column defining the factor (e.g. "Treatment") in code{phenodata} is not 
  an R factor, use \code{pData(targets)$Treatment = factor(pData(targets)$Treatment)}
  to convert the character values into an actual R factor with appropriate levels.
}
\seealso{
  Methods
  \code{\link[GOexpress:subset_scores]{GOexpress::subset_scores}},
  \code{\link[biomaRt:getBM]{biomaRt::getBM}},
  \code{\link[randomForest:randomForest]{randomForest::randomForest}}, 
  \code{\link[stats:oneway.test]{stats::oneway.test}},
  and
  \href{https://github.com/kevinrue/GOexpress/wiki/Typical-analysis-pipeline}{tutorial} on GitHub.
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
raw_result = GOexpress::GO_analyse(expr_data=log2cpm, phenodata=targets, f="Treatment")
}

\keyword{ GOexpress }
\keyword{ randomForest }
\keyword{ anova }
\keyword{ gene }
\keyword{ expression }
\keyword{ clustering }
\keyword{ ontology }
