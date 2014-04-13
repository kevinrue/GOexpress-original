\name{subset.scores}
\alias{subset.scores}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
subset.scores(result, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{result}{
%%     ~~Describe \code{result} here~~
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
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
function (result, ...) 
{
    print(list(...))
    print(names(list(...)))
    filters = list(...)
    filtered = data.frame(row.names = result$scores$go_id)
    for (filter in names(list(...))) {
        if (!filter \%in\% colnames(result$scores)) {
            stop(filter, " is not a valid column name.")
        }
        if (filter \%in\% c("gene_count", "gene_count", "ave.F.score")) {
            cat(filter, " equal or more than:", filters[[filter]], 
                fill = TRUE)
            filtered[, filter] = result$scores[, filter] >= filters[filter]
        }
    }
    filtered$merge = apply(X = filtered, MARGIN = 1, FUN = all)
    list(scores = subset(result$score, filtered$merge), mapping = result$mapping)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
