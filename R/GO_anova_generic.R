# This script will take
## An expression dataset (using ensembl_id or probeset identifiers)
## A table of phenotypic data of samples
## The phenotypic factor expected to discriminate samples

# The script will evaluate the power of each expressed gene to 
# discriminate the samples according to the expected phenotypic factor

# The script will summarise genes annotated to a same GO term (Biological Process)
# to evaluate the power of each GO term to discriminate the samples according 
# to the expected phenotypic factor

# The script will return
## A table ranking GO terms according to their power to discriminate samples
## A table mapping  GO terms and genes to one another
## A table ranking genes according to their power to discriminate samples
## The factor used to calculate scores

# Downstream analyses
## Function to filter the resulting scores to GO terms passing criteria
## Function to list the genes annotated to a GO term
## Function to see the scores of the genes annotated to a GO term
## Function to cluster samples according to a list of genes
## Function to heatmap the samples according to a list of genes
## Function to see the distribution of scores for all GO terms (histogram)
## Function to see the quantiles of scores across all GO terms
## Function to see quantiles of scores for all GO terms (how many GO terms have score > S)
## Function to filter the results for GO terms with N genes or more
## Two function to visualise the expression profile of a gene across a X-variable and 
### grouped by a Y-factor (plot by gene_id or by gene_symbol)
## Function to visualise the univariate effect of each factor on the grouping of samples


# Dependencies:
## Internet connection (biomaRt)
## library(biomaRt) # gene annotation
## library(stringr) # pattern extraction


GO_anova = function(expr_data, phenodata, f, biomart_dataset="", microarray="",
                    adj.P.method = "BH", FDR=0.05){
  # if the user did not give a dataset name
  if (biomart_dataset == ""){
    # if the user did not give a microarray value
    if (microarray == ""){
      # automatically detect both
      # fetch the first gene id in the given expression dataset
      sample_gene = rownames(expr_data)[1]
      # if the gene id has not an identifiable ensembl gene id prefix
      mart = mart_from_ensembl(sample_gene)
      if (!class(mart) == "Mart"){
        # if the gene id has an identifiable microarray gene id prefix
        microarray_match = microarray_from_probeset(sample_gene)
        if (!is.null(nrow(microarray_match))){
          # connect to biomart and set the microarray variable
          cat("Looks like microarray data.", fill=TRUE)
          cat("Loading dataset", microarray_match$dataset, "for microarray",
              microarray_match$microarray, "...", fill=TRUE)
          microarray = microarray_match$microarray
          biomart_dataset = microarray_match$dataset
          mart = useMart(biomart="ensembl", dataset=biomart_dataset)
        }
        # if the gene id does not have an identifiable microarray gene id prefix
        else{
          # stop the program and throw an error
          stop("Cannot guess origin of dataset:
Please use \"biomart_dataset=\" and/or \"microarray=\" arguments.")
        }
      }
      # if the gene id has an identifiable ensembl gene id prefix
      # then, the connection to the mart is already established
      # leave microarray to "" to imply that we don't work with microarray data
    }
    # if the user gave a microarray name
    else{
      # check if it exists
      if (!microarray %in% microarray2dataset$microarray){
        stop("Invalid microarray value. See data(microarray2dataset)")
      }
      # if it is unique to a dataset (some microarray have the same column name
      if(sum(microarray2dataset$microarray == microarray) == 1){
        biomart_dataset = microarray2dataset[microarray2dataset$microarray == microarray, "dataset"]
        cat("Loading requested microarray", microarray, "from detected dataset", biomart_dataset, "...", fill=TRUE)
        mart = useMart(biomart="ensembl", dataset=biomart_dataset)
        # Leave microarray to the current valid value
      }
      # if the microarray does not exist in the dataset
      else if(sum(microarray2dataset$microarray == microarray) == 0){
        stop("Microarray name not recognised. See data(microarray2dataset).")
      }
      # if the microarray name exists in multiple datasets
      else{
        cat("Multiple datasets possible:", fill=TRUE)
        print(microarray2dataset[microarray2dataset$microarray == microarray, c("dataset", "microarray")])
        stop("Cannot guess dataset.
Please use \"biomart_dataset=\" argument.")
      }  
    }
  }
  # if the user gave a biomart_dataset value
  else{
    # Check that it exists
    if (!biomart_dataset %in% prefix2dataset$dataset){
      stop("Invalid biomart_dataset value. See data(prefix2dataset)")
    }
    cat("Using biomart dataset", biomart_dataset, fill=TRUE)
    # if the user did not give a microarray name
    if (microarray == ""){
      # Check if looks like microarray
      # fetch the first gene id in the given expression dataset
      sample_gene = rownames(expr_data)[1]
      microarray_match = microarray_from_probeset(sample_gene)
      # if the data matches a known microarray pattern
      if (!is.null(nrow(microarray_match))){
        # connect to biomart and set the microarray variable
        cat("Looks like microarray data.", fill=TRUE)
        # if the dataset/microarray pair exists
        if (microarray_match$dataset == biomart_dataset){
          cat("Loading annotations for microarray",
              microarray_match$microarray, "...", fill=TRUE)
          mart = useMart(biomart="ensembl", dataset=biomart_dataset)
          microarray = microarray_match$microarray
          print(mart)
        }
        # if the dataset/microarray pair does no exist
        else{
          # The dataset exists, the data matches a microarray
          # but not a microarray of the dataset
          cat("Detected microarray", microarray_match$microarray,
               "inexisting in requested dataset", biomart_dataset,
               ". Possible datasets are:")
          return(microarray_match)
        }
      }
      # if the data does not match a microarray pattern
      else{
        cat(sample_gene, "gene identifier in expression data cannot\
  be resolved to a microarray. Assuming ensembl gene identifiers.", fill=TRUE)
      }
      # If it does not look like microarray
      # assume it is ensembl annotations
      # therefore do nothing more
      # in both cases load the requested mart dataset
      cat("Loading requested dataset", biomart_dataset, "...", fill=TRUE)
      mart = useMart(biomart="ensembl", dataset=biomart_dataset)
    }
    # if the user gave a microarray name
    else{
      # Check that the pair dataset/microarray exists
      if (!biomart_dataset %in% microarray2dataset[microarray2dataset$microarray == microarray, "dataset"]){
        stop("There is no microarray ", microarray, " in dataset ", biomart_dataset)
      }
      cat("Loading requested microarray", microarray, "from biomart dataset", biomart_dataset, fill=TRUE)
      mart = useMart(biomart="ensembl", dataset=biomart_dataset)
    }
  }
  print(mart)
  #  if working with ensembl gene identifiers
  if (microarray == ""){
    # Prepare a mapping table between gene identifiers and GO terms
    cat("Fetching ensembl_gene_id/GO_id mappings from BioMart ...", fill=TRUE)
    GO_genes = getBM(attributes=c("ensembl_gene_id", "go_id"), mart=mart)
  }
  # if working with microarray probesets
  else{
    # Prepare a mapping table between gene identifiers and GO terms
    cat("Fetching probeset/GO_id mappings from BioMart ...", fill=TRUE)
    GO_genes = getBM(attributes=c(microarray, "go_id"), mart=mart)
  }
  # Rename the first column which could be ensembl_id or probeset_id
  colnames(GO_genes)[1] = "gene_id"
  # Remove over 1,000 rows where the go_id is ""
  GO_genes = GO_genes[GO_genes$go_id != "",]
  # Remove rows where the gene_id is "" (happens)
  GO_genes = GO_genes[GO_genes$gene_id != "",]
  # Prepare a table of all the GO termsin BioMart (even if no gene is annotated to it)
  cat("Fetching GO_terms description from BioMart ...", fill=TRUE)
  all_GO = getBM(attributes=c("go_id", "name_1006", "namespace_1003"),
                 mart=mart)
  # Remove the GO terms which is ""
  all_GO = subset(x=all_GO, subset=go_id != "")
  # Calculate the F.value and p.value of ANOVA for each ensembl id in the expression dataset
  cat("Calculating one-way ANOVA on factor", f,"for", nrow(expr_data), "genes. This may take a few minutes ... (about 2min for 12,000 genes)", fill=TRUE)
  res_anova = data.frame("F.value"= apply(X=expr_data,
                                          MARGIN=1,
                                          FUN=function(x){oneway.test(formula=expr~group,
                                                                      data=cbind(expr=x, group=pData(phenodata)[,f]))$statistic}),
                         "p.value"=apply(X=expr_data,
                                         MARGIN=1,
                                         FUN=function(x){oneway.test(formula=expr~group,
                                                                     data=cbind(expr=x, group=pData(phenodata)[,f]))$p.value}))
  # Correct for multiple testing (default: BH)
  res_anova$FDR = p.adjust(p=res_anova$p.value, method=adj.P.method)
  # Set to 0 the F.value of non-significant genes
  res_anova$F.value[res_anova$FDR > FDR] = 0
  # Summary statistics by GO term
  ## Merge the table mapping GOterm to genes with the ANOVA results of each gene (twice)
  # First merge the tables while keeping all gene/GO mappings, even for genes absent of the dataset
  # This will allow average F values to be calculated on the basis of all ensembl genes annotated to 
  # the GO term, even if not in the dataset (genes absent are considered non-significant)
  GO_gene_anova_all = merge(x=GO_genes, y=res_anova, by.x="gene_id", by.y="row.names", all.x=TRUE)
  # The merge will leave genes absent from the dataset with NA as F.value and p.value. 
  # For the statistics to work, we only need to replace F.value by 0 where there are NAs
  GO_gene_anova_all[is.na(GO_gene_anova_all$F.value),]$F.value = 0
  # Second, merge the tables keeping only the genes present in the dataset, this value is not used elsewhere
  # in the program, but can be valued by the user
  GO_gene_anova_data = merge(x=GO_genes, y=res_anova, by.x="gene_id", by.y="row.names")
  # Results can now be summarised by aggregating rows with same GOterm
  # Appends gene annotations to rows of res_anova
  cat("Fetching gene description from BioMart ...", fill=TRUE)
  #  if working with ensembl gene identifiers
  if (microarray == ""){
    genes_anova = merge(x=res_anova, all.x=TRUE,
                        y=getBM(attributes=c("ensembl_gene_id", "external_gene_id", "description"),
                                filters="ensembl_gene_id",
                                values=rownames(res_anova),
                                mart=mart),
                        by.x="row.names",
                        by.y="ensembl_gene_id")
  }
  # if working with microarray probesets
  else{
    # Prepare a mapping table between gene identifiers and GO terms
    genes_anova = merge(x=res_anova, all.x=TRUE, 
                        y=getBM(attributes=c(microarray, "external_gene_id", "description"),
                                filters=microarray,
                                values=rownames(res_anova),
                                mart=mart),
                        by.x="row.names",
                        by.y=microarray)
    # In the case of microarray, probesets can be annotated to multiple gene symbols
    # for each probeset, keep only the first row (we're losing one possible annotation here)
    genes_anova = genes_anova[ !duplicated(genes_anova$Row.names), ]
  }
  # Rank the genes by increasing FDR (should match increasing F.value)
  genes_anova = genes_anova[order(genes_anova$FDR),]
  # Put the ensembl identifier back as the row name
  rownames(genes_anova) = genes_anova$Row.names
  genes_anova$Row.names = NULL
  cat("Merging score into result table ...", fill=TRUE)
  # Number of significant genes in dataset to each GO term
  GO_scores = merge(x=aggregate(F.value~go_id, data=GO_gene_anova_data, FUN=function(x){sum(x != 0)}), y=all_GO, by="go_id", all.y=TRUE)
  colnames(GO_scores)[2] = "sig_count"
  # Total number of genes in the dataset annotated to each GO term
  GO_scores = merge(x=aggregate(gene_id~go_id, data=GO_gene_anova_data, FUN=length), y=GO_scores, by="go_id", all.y=TRUE)
  colnames(GO_scores)[2] = "data_count"
  # Total number of genes annotated to each GO term in BioMart (not necessarily in dataset)
  GO_scores = merge(x=aggregate(F.value~go_id, data=GO_gene_anova_all, FUN=length), y=GO_scores, by="go_id", all.y=TRUE)
  colnames(GO_scores)[2] = "total_count"
  ## Average F value (denominator being the total of genes by GO term in BioMart) being tested
  # (+) robust for GO terms with several genes (5 minimum advised, 10 was found robust, gene counts per GO term below)
  GO_scores = merge(x=aggregate(F.value~go_id, data=GO_gene_anova_all, FUN=mean), y=GO_scores, by="go_id", all.y=TRUE)
  colnames(GO_scores)[2] = "ave.F.score"  
  # Notes of other metrics tested:
  ## Sum.F.values: (-) biased toward general GO terms annotated for many thousands of genes (e.g. "protein binding")
  ## Max.F.values: (+) insensitive to number of genes annotated for GO term
  #                (-) many GO terms sharing the same gene are tied (-) not a robust metric of GO term
  # Most top ranked GO terms according to the average F value contain a single gene
  # But this bias can easily be attenuated by filtering for GO terms with a minimal number of genes
  # Rank the GO terms by decreasing average F value
  GO_scores = GO_scores[order(GO_scores$ave.F.score, decreasing=TRUE),]
  # Return the results of the analysis
  return(list(scores=GO_scores, mapping=GO_genes, anova=genes_anova, factor=f))
}

mart_from_ensembl = function(sample_gene){
  # If the gene id starts by "ENS" (most cases, except 3 handled separately below)
  if (length(grep(pattern="^ENS", x=sample_gene))){
    # Extract the full prefix
    prefix = str_extract(sample_gene, "ENS[[:upper:]]+")
    # If the ENS* prefix is in the table 
    if (prefix %in% prefix2dataset$prefix){
      # load the corresponding biomart dataset
      cat("Looks like ensembl gene identifier", fill=TRUE)
      cat("Loading dataset", prefix2dataset[prefix2dataset$prefix == prefix,]$dataset,
          "...", fill=TRUE)
      return(useMart(biomart="ensembl",
                     dataset=prefix2dataset[prefix2dataset$prefix == prefix,]$dataset))
    }
    # Otherwise return FALSE
    else{
      cat("Did not recognise a valid ensembl gene identifier.", fill=TRUE)
      return(FALSE)
    }
  }
  # If the gene id starts with "WBgene"
  else if (length(grep(pattern="^WBGene", x=sample_gene))) {
    # load the corresponding biomart dataset
    cat("Looks like ensembl gene identifier", fill=TRUE)
    cat("Loading dataset celegans_gene_ensembl ...", fill=TRUE)
    return(useMart(biomart="ensembl", dataset="celegans_gene_ensembl"))
  }
  # If the gene id starts with "FBgn"
  else if (length(grep(pattern="^FBgn", x=sample_gene))) {
    # load the corresponding biomart dataset
    cat("Looks like ensembl gene identifier", fill=TRUE)
    cat("Loading dataset dmelanogaster_gene_ensembl ...", fill=TRUE)
    return(useMart(biomart="ensembl", dataset="dmelanogaster_gene_ensembl"))
  }
  # If the gene id starts with "Y"
  else if (length(grep(pattern="^Y", x=sample_gene))) {
    # load the corresponding biomart dataset
    cat("Looks like ensembl gene identifier", fill=TRUE)
    cat("Loading dataset scerevisiae_gene_ensembl ...", fill=TRUE)
    return(useMart(biomart="ensembl", dataset="scerevisiae_gene_ensembl"))
  }
  # If the gene id does not match any known ensembl gene id prefix, return an error and stop
  else{
    cat("Did not recognise an ensembl gene identifier.", fill=TRUE)
    return(FALSE)
  }
}

microarray_from_probeset = function(sample_gene){
  matches = c()
  # For each pattern thought to be unique to a microarray
  for (pattern in microarray2dataset$prefix[microarray2dataset$unique]){
    # if the pattern matches the sample gene
    if (length(grep(pattern=pattern, x=sample_gene))){
      # add the pattern to a vector 
      matches = c(matches, pattern)
    }
  }
  # if the vector length is at least 1 (unlikely to ever be more than 1 if patterns do not overlap)
  if (length(matches)){
    # return (dataset, microarray) to the main function to build mapping tables
    return(microarray2dataset[microarray2dataset$prefix == matches[1],
                                           c("dataset","microarray")])
  }
  # If the sample gene was not recognised in the unique ones,
  # check whether it may be an ambiguous identifier
  # For each unique pattern known to be found in multiple microarrays
  for (pattern in unique(microarray2dataset$prefix[!microarray2dataset$unique])){
    # if the pattern matches the sample gene
    if (length(grep(pattern=pattern, x=sample_gene))){
      # add the pattern to a vector 
      matches = c(matches, pattern)
    }
  }
  # if vector contains at least 1 pattern
  if (length(matches)){
    # print sample, first pattern, and list of possible microarray
    cat(sample_gene, "matches pattern", matches[1], "found in multiple microarrays:", fill=TRUE)
    print(microarray2dataset[microarray2dataset$prefix == matches[1], 
                              c("dataset","microarray")], row.names=FALSE)
    return(FALSE)
  }
  # if no known microarray pattern matches, return FALSE
  else{
    cat("Did not recognise microarray data.", fill=TRUE)
    return(FALSE)
  }
}
