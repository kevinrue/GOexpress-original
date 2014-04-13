\name{GO_anova}
\alias{GO_anova}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
GO_anova(expr_data, phenodata, f, biomart_dataset = "", adj.P.method = "BH")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{expr_data}{
%%     ~~Describe \code{expr_data} here~~
}
  \item{phenodata}{
%%     ~~Describe \code{phenodata} here~~
}
  \item{f}{
%%     ~~Describe \code{f} here~~
}
  \item{biomart_dataset}{
%%     ~~Describe \code{biomart_dataset} here~~
}
  \item{adj.P.method}{
%%     ~~Describe \code{adj.P.method} here~~
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
function (expr_data, phenodata, f, biomart_dataset = "", adj.P.method = "BH") 
{
    library(biomaRt, quietly = TRUE)
    cat("Connecting to appropriate biomaRt dataset...", fill = TRUE)
    mart = getMart(expr_data, biomart_dataset)
    print(mart)
    cat("Fetching ensembl gene/GO term mapping from biomaRt...", 
        fill = TRUE)
    GO_genes = getBM(attributes = c("ensembl_gene_id", "go_id"), 
        filters = "ensembl_gene_id", values = rownames(expr_data), 
        mart = mart, uniqueRows = TRUE)
    GO_genes = subset(GO_genes, go_id != "")
    cat("Calculating ANOVA for", nrow(expr_data), "genes. This may take a few minutes... (about 2min for 12,000 genes)", 
        fill = TRUE)
    res_anova = data.frame(F.value = apply(X = log2_cpm_norm_no0H, 
        MARGIN = 1, FUN = function(x) {
            oneway.test(formula = expr ~ group, data = cbind(expr = x, 
                group = pData(targets_no0H)[, f]))$statistic
        }), p.value = apply(X = log2_cpm_norm_no0H, MARGIN = 1, 
        FUN = function(x) {
            oneway.test(formula = expr ~ group, data = cbind(expr = x, 
                group = pData(targets_no0H)[, f]))$p.value
        }))
    res_anova$FDR = p.adjust(p = res_anova$p.value, method = adj.P.method)
    res_anova$F.value[res_anova$FDR > 0.05] = 0
    GO_gene_anova = merge(x = GO_genes, y = res_anova, by.x = "ensembl_gene_id", 
        by.y = "row.names")
    cat("Fetching GO terms description from biomaRt...", fill = TRUE)
    all_GO = getBM(attributes = c("go_id", "name_1006", "namespace_1003"), 
        filters = "ensembl_gene_id", values = rownames(log2_cpm_norm_no0H), 
        mart = mart, uniqueRows = TRUE)
    all_GO = subset(all_GO, go_id != "")
    cat("Merging score into result table...", fill = TRUE)
    GO_scores = merge(x = aggregate(F.value ~ go_id, data = GO_gene_anova, 
        FUN = mean), y = all_GO, by = "go_id")
    colnames(GO_scores)[2] = "ave.F.score"
    GO_scores = merge(x = aggregate(F.value ~ go_id, data = GO_gene_anova, 
        FUN = function(x) {
            sum(x != 0)
        }), y = GO_scores, by = "go_id")
    colnames(GO_scores)[2] = "sig_count"
    GO_scores = merge(x = aggregate(ensembl_gene_id ~ go_id, 
        data = GO_gene_anova, FUN = length), y = GO_scores, by = "go_id")
    colnames(GO_scores)[2] = "gene_count"
    GO_scores = GO_scores[order(GO_scores$ave.F.score, decreasing = TRUE), 
        ]
    return(list(scores = GO_scores, mapping = GO_genes))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
