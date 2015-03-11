\name{subset_scores}

\alias{subset_scores}

\title{
Returns a filtered list from GO_analyse results.
}

\description{
Given an output variable from a \code{GO_analyse} analysis and a set of valid
filters and thresholds, returns an identically formatted list keeping only the
rows of the score table passing all the filters.
}

\usage{
subset_scores(result, ...)

# Suggested use:
# subset_scores(result, total_count=5, namespace_1003="biological_process")
# synonym to 
# subset_scores(result, total=5, namespace="BP")
}

\arguments{
    \item{result}{
    The output of the \code{GO_analyse()} function.
    }
    \item{\dots}{
    Additional pairs of filter and threshold values in the format
    "filter=threshold". Filters must be valid names from
    colnames(result$GO).
    }
}

\details{
It is critical to filter out GO terms with very few genes (e.g. less than 5
genes), as the scoring function is biased for those GO terms (see UsersGuide).

In addition, it is useful to retain only GO terms of one type (i.e. namespace)
among the three possible: "biological process", "molecular function", and
"cellular component".

Suggested filters:
\tabular{ll}{
total_count, total: \tab Filter keeping only GO terms associated with at
least the given count of genes in the annotations provided.\cr
namespace, namespace_1003: \tab Filter keeping only the GO terms of a given
type. Valid values are "biological_process", "molecular_function", and
"cellular_component". Abbreviations "BP", "MF", and "CC" are also accepted.\cr
p.val: \tab Filter keeping only the GO terms with P-value lower or equal to
the given cutoff (This filter is only applicable after the use of the
\code{pValue_GO} function).\cr
}

Other filters:
\tabular{ll}{
data_count, data: \tab Filter keeping only GO terms associated with at
least the given count of genes, present in the \code{ExpressionSet}.\cr

ave_rank: \tab Filter keeping only GO terms with an ave_rank value equal
or lower than the given cutoff (average of the rank of all genes annotated
to the GO term).\cr

ave_score: \tab Filter keeping only GO terms with an ave_score value equal
or higher than the given cutoff (average of the score of all genes annotated
to the GO term) Scores are the mean decrease in Gini index for the random
forest, or the F-value for the ANOVA approach.\cr
}

}

\value{
A list formatted identically to the results of the analysis, but restricted to
the gene ontologies passing the given filters, and the genes mapped to those
remaining ontologies (i.e. all genes not associated with an ontology are
discarded in the output object). Additionally, only mapping between the
remaining gene ontologies will be retained in the output object, 
to limit memory usage.
}

\author{
    Kevin Rue-Albrecht
}

\seealso{
    Method \code{\link[GOexpress:GO_analyse]{GO_analyse}}.
}

\examples{
# load the sample output data
data(AlvMac_results.pVal)

# have an overview of the result variable
str(AlvMac_results.pVal)

# filter for Biological Processes associated with 5+ genes and <=0.05 P-value
filtered_results <- subset_scores(
    result=AlvMac_results.pVal, total_count=5, p.val=0.05,
    namespace="BP")

# have an overview of the filtered result variable
str(filtered_results)
}

\keyword{ GOexpress }
