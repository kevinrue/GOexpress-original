# This script will take
## An expression dataset (requires ensembl-compatible identifiers)
## A table of phenotypic data of samples
## The phenotypic factor expected to discriminate samples

# The script will evaluate the power of each expressed gene to 
# discriminate the samples according to the expected phenotypic factor

# The script will summarise genes annotated to a same GO term (Biological Process)
# to evaluate the pwoer of each GO term to discriminate the samples according 
# to the expected phenotypic factor

# The script will return
## A table ranking GO terms according to their power to discriminate samples
## A table mapping  GO terms and genes to one another
## A table describing all GO terms annotated for genes in the dataset


# Downstream analyses
## Contains a function to list the genes annotated to a GO term
## Contains a function to cluster samples according to a list of genes
## Contains a function to heatmap the samples according to a list of genes
## Contains two functions to see the distribution of scores for all GO terms (density + histogram)
## Contains a function to see quantiles of scores for all GO terms (how many GO terms have score > S)
## Contains a function to filter the results for GO terms with N genes or more

# Dependencies:
## Internet connection (biomaRt)
## library(biomaRt) # gene annotation
## library(stringr) # pattern extraction


GO_anova = function(expr_data, phenodata, f, biomart_dataset="", adj.P.method = "BH"){
  library(biomaRt, quietly=TRUE)
  # Connect to the appropriate biomaRt for the expression data or the user specified dataset
  cat("Connecting to appropriate biomaRt dataset...", fill=TRUE)
  mart = getMart(expr_data, biomart_dataset)
  # Information for the user
  print(mart)
  # Prepare a table mapping the ensembl ids in the expression dataset to GO terms (biological processes)
  cat("Fetching ensembl gene/GO term mapping from biomaRt...", fill=TRUE)
  GO_genes = getBM(attributes=c("ensembl_gene_id", "go_id"),
                   filters="ensembl_gene_id",
                   values=rownames(expr_data),
                   mart=mart,
                   uniqueRows=TRUE)
  # Remove over 1,000 rows where the go_id is ""
  GO_genes = subset(GO_genes, go_id != "")
  # Calculate the F.value and p.value of ANOVA for each ensembl id in the expression dataset
  cat("Calculating ANOVA for", nrow(expr_data), "genes. This may take a few minutes... (about 2min for 12,000 genes)", fill=TRUE)
  res_anova = data.frame("F.value"= apply(X=log2_cpm_norm_no0H,
                                          MARGIN=1,
                                          FUN=function(x){oneway.test(formula=expr~group,
                                                                      data=cbind(expr=x, group=pData(targets_no0H)[,f]))$statistic}),
                         "p.value"=apply(X=log2_cpm_norm_no0H,
                                         MARGIN=1,
                                         FUN=function(x){oneway.test(formula=expr~group,
                                                                     data=cbind(expr=x, group=pData(targets_no0H)[,f]))$p.value}))
  # Correct for multiple testing (default: BH)
  res_anova$FDR = p.adjust(p=res_anova$p.value, method=adj.P.method)
  # Set to 0 the F.value of non-significant genes
  res_anova$F.value[res_anova$FDR > 0.05] = 0
  # Summary statistics by GO term
  ## Merge the table mapping GOterm to genes with the ANOVA results of each gene
  GO_gene_anova = merge(x=GO_genes, y=res_anova, by.x="ensembl_gene_id", by.y="row.names")
  # Results can now be summarised by aggregating rows with same GOterm
  ## Prepare a table of all the GO terms linked to at least a gene in the expression data
  cat("Fetching GO terms description from biomaRt...", fill=TRUE)
  all_GO = getBM(attributes=c("go_id", "name_1006", "namespace_1003"),
                 filters="ensembl_gene_id",
                 values=rownames(log2_cpm_norm_no0H),
                 mart=mart,
                 uniqueRows=TRUE)
  # Remove one of the GO terms which is ""
  all_GO = subset(all_GO, go_id != "")
  cat("Merging score into result table...", fill=TRUE)
  ## Average F value (+) robust for GO terms with several genes (5 minimum advised, 10 was found robust, gene counts per GO term below)
  GO_scores = merge(x=aggregate(F.value~go_id, data=GO_gene_anova, FUN=mean), y=all_GO, by="go_id")
  colnames(GO_scores)[2] = "ave.F.score"
  # Notes of other metrics tested:
  ## Sum.F.values: (-) biased toward general GO terms annotated for many thousands of genes (e.g. "protein binding")
  ## Max.F.values: (+) insensitive to number of genes annotated for GO term (-) many GO terms sharing the same gene are tied (-) not a robust metric of GO term
  # Most top ranked GO terms contain a single gene
  # The 
  # Number of significant genes by GO term
  GO_scores = merge(x=aggregate(F.value~go_id, data=GO_gene_anova, FUN=function(x){sum(x != 0)}), y=GO_scores, by="go_id")
  colnames(GO_scores)[2] = "sig_count"
  # Total number of genes annotated by GO term
  GO_scores = merge(x=aggregate(ensembl_gene_id~go_id, data=GO_gene_anova, FUN=length), y=GO_scores, by="go_id")
  colnames(GO_scores)[2] = "gene_count"
  # Rank the GO terms by decreasing average F value
  GO_scores = GO_scores[order(GO_scores$ave.F.score, decreasing=TRUE),]
  # Return the results of the analysis
  return(list(scores=GO_scores, mapping=GO_genes))
}

hist.scores = function(result){
  hist(result$scores$ave.F.score)
}

quantiles.scores = function(result, quartiles=c(FALSE, TRUE)){
  # If user changes to "quartiles=TRUE" then the following will be true
  if (quartiles[1]){
    quantile(x=result$scores$ave.F.score)
  }
  # Default is to give range of top 10%, 5%, 1%, 0.1% and 0.01%
  else{
    quantile(x=result$scores$ave.F.score, c(0.9, 0.95, 0.99, 0.999, 0.9999))
  }
}

show.scores = function(result){
  result$scores
}

subset.scores = function(result, ...){
  print(list(...))
  print(names(list(...)))
  #print(head(result$scores))
  # Save the list of filter and value for easier referencing
  filters = list(...)
  # prepares a table where the filters will be saved
  filtered = data.frame(row.names=result$scores$go_id)
  # For each filter
  for (filter in names(list(...))){
    # Check that it is a valid filter
    if (!filter %in% colnames(result$scores)){
      stop(filter, " is not a valid column name.")
    }
    # Save the filter status of each row for this filter
    ## Filters where the value should be higher than the given threshold
    if (filter %in% c("gene_count", "gene_count", "ave.F.score")){
      cat(filter, " equal or more than:", filters[[filter]], fill=TRUE)
      filtered[,filter] = result$scores[,filter] >= filters[filter]
    }
  }
  #
  filtered$merge = apply(X=filtered, MARGIN=1, FUN=all)
  #
  list(scores=subset(result$score, filtered$merge), mapping=result$mapping)
}

getMart = function(expr_data, biomart_dataset){
  # This function returns a mart from biomaRt
  ## either the dataset specified by the user if any
  ## or the dataset corresponding to the species guessed from an gene id in the expression data
  # If the biomart dataset is specified by user
  if (biomart_dataset != ""){
    # load the biomart dataset requested by the user (let it crash if incorrect dataset specified)
    return(useMart(biomart="ensembl", biomart_dataset))
  }
  # Otherwise, guess species from the expression data
  else {
    # Take the first gene id
    sample_gene = rownames(expr_data)[1]
    # If the gene id starts by "ENS" (most cases)
    if (length(grep(pattern="^ENS", x=sample_gene))){
      library(stringr)
      # Extract the full prefix
      prefix = str_extract(sample_gene, "ENS[[:upper:]]+")
      # load the precomputed table mapping prefix to biomart database
      load("prefix2dataset.RData")
      # If the prefix is in the table 
      if (prefix %in% prefix2dataset$prefix){
        # load the corresponding biomart dataset
        return(useMart(biomart="ensembl", dataset=prefix2dataset[prefix2dataset$prefix == prefix,]$dataset))
      }
      # Otherwise return an error and stop
      else{
        stop(prefix, " is not recognised as an ensembl gene id prefix.")
      }
    }
    # If the gene id starts with "WBgene"
    else if (length(grep(pattern="^WBGene", x=sample_gene))) {
      # load the corresponding biomart dataset
      return(useMart(biomart="ensembl", dataset="celegans_gene_ensembl"))
    }
    # If the gene id starts with "FBgn"
    else if (length(grep(pattern="^FBgn", x=sample_gene))) {
      # load the corresponding biomart dataset
      return(useMart(biomart="ensembl", dataset="dmelanogaster_gene_ensembl"))
    }
    # If the gene id starts with "Y"
    else if (length(grep(pattern="^Y", x=sample_gene))) {
      # load the corresponding biomart dataset
      retur(useMart(biomart="ensembl", dataset="scerevisiae_gene_ensembl"))
    }
    # If the gene id does not match any known ensembl gene id prefix, return an error and stop
    else{
      stop(sample_gene, " is not recognised as an ensembl gene id.")
    }
  }
}