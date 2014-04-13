\name{getMart}
\alias{getMart}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
getMart(expr_data, biomart_dataset)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{expr_data}{
%%     ~~Describe \code{expr_data} here~~
}
  \item{biomart_dataset}{
%%     ~~Describe \code{biomart_dataset} here~~
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
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (expr_data, biomart_dataset) 
{
    if (biomart_dataset != "") {
        return(useMart(biomart = "ensembl", biomart_dataset))
    }
    else {
        sample_gene = rownames(expr_data)[1]
        if (length(grep(pattern = "^ENS", x = sample_gene))) {
            library(stringr)
            prefix = str_extract(sample_gene, "ENS[[:upper:]]+")
            load("prefix2dataset.RData")
            if (prefix \%in\% prefix2dataset$prefix) {
                return(useMart(biomart = "ensembl", dataset = prefix2dataset[prefix2dataset$prefix == 
                  prefix, ]$dataset))
            }
            else {
                stop(prefix, " is not recognised as an ensembl gene id prefix.")
            }
        }
        else if (length(grep(pattern = "^WBGene", x = sample_gene))) {
            return(useMart(biomart = "ensembl", dataset = "celegans_gene_ensembl"))
        }
        else if (length(grep(pattern = "^FBgn", x = sample_gene))) {
            return(useMart(biomart = "ensembl", dataset = "dmelanogaster_gene_ensembl"))
        }
        else if (length(grep(pattern = "^Y", x = sample_gene))) {
            retur(useMart(biomart = "ensembl", dataset = "scerevisiae_gene_ensembl"))
        }
        else {
            stop(sample_gene, " is not recognised as an ensembl gene id.")
        }
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
