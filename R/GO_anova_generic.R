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


GO_anova = function(expr_data, phenodata, f, biomart_dataset="", adj.P.method = "BH", FDR=0.05){
  # Connect to the appropriate biomaRt for the expression data or the user specified dataset
  mart = get_mart_dataset(expr_data, biomart_dataset)
  # Prepare a table mapping all ensembl ids in the BioMart dataset to GO terms (biological processes)
  # Not just the ensembl ids in the data, but all in BioMart (will affect the GO term average score)
  cat("Fetching ensembl_gene/GO_term mapping from BioMart ...", fill=TRUE)
  GO_genes = getBM(attributes=c("ensembl_gene_id", "go_id"),
                   mart=mart)
  # Remove over 1,000 rows where the go_id is ""
  GO_genes = GO_genes[GO_genes$go_id != "",]
  # Prepare a table of all the GO termsin BioMart (even if no gene is annotated to it)
  cat("Fetching GO_terms description from BioMart ...", fill=TRUE)
  all_GO = getBM(attributes=c("go_id", "name_1006", "namespace_1003"),
                 mart=mart)
  # Remove the GO terms which is ""
  all_GO = subset(x=all_GO, subset=go_id != "")
  # Calculate the F.value and p.value of ANOVA for each ensembl id in the expression dataset
  cat("Calculating one-way ANOVA on factor \"", f,"\" for", nrow(expr_data), "genes. This may take a few minutes ... (about 2min for 12,000 genes)", fill=TRUE)
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
  GO_gene_anova_all = merge(x=GO_genes, y=res_anova, by.x="ensembl_gene_id", by.y="row.names", all.x=TRUE)
  # The merge will leave genes absent from the dataset with NA as F.value and p.value. 
  # For the statistics to work, we only need to replace F.value by 0 where there are NAs
  GO_gene_anova_all[is.na(GO_gene_anova_all$F.value),]$F.value = 0
  # Second, merge the tables keeping only the genes present in the dataset, this value is not used elsewhere
  # in the program, but can be valued by the user
  GO_gene_anova_data = merge(x=GO_genes, y=res_anova, by.x="ensembl_gene_id", by.y="row.names")
  # Results can now be summarised by aggregating rows with same GOterm
  # Appends gene annotations to rows of res_anova
  cat("Fetching gene description from BioMart ...", fill=TRUE)
  genes_anova = merge(x=res_anova,
                      y=getBM(attributes=c("ensembl_gene_id", "external_gene_id", "description"),
                              filters="ensembl_gene_id",
                              values=rownames(res_anova),
                              mart=mart),
                      by.x="row.names",
                      by.y="ensembl_gene_id")
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
  GO_scores = merge(x=aggregate(ensembl_gene_id~go_id, data=GO_gene_anova_data, FUN=length), y=GO_scores, by="go_id", all.y=TRUE)
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

get_mart_dataset = function(expr_data, biomart_dataset){
  # This function returns a mart from biomaRt
  ## either the dataset specified by the user if any
  ## or the dataset corresponding to the species guessed from an gene id in the expression data
  # If the biomart dataset is specified by user
  if (biomart_dataset != ""){
    cat("Connecting to requested ensembl dataset", biomart_dataset, "...", fill=TRUE)
    # load the biomart dataset requested by the user (let it crash if incorrect dataset specified)
    return(useMart(biomart="ensembl", biomart_dataset))
  }
  # Otherwise, guess species from the expression data
  else {
    cat("Guessing appropriate BioMart dataset ...", fill=TRUE)
    # Take the first gene id
    sample_gene = rownames(expr_data)[1]
    cat("First gene identifier in expression dataset:", sample_gene, fill=TRUE)
    # If the gene id starts by "ENS" (most cases)
    if (length(grep(pattern="^ENS", x=sample_gene))){
      # Extract the full prefix
      prefix = str_extract(sample_gene, "ENS[[:upper:]]+")
      # If the prefix is in the table 
      if (prefix %in% prefix2dataset$prefix){
        cat("Looks like you need ensembl dataset:", 
            prefix2dataset[prefix2dataset$prefix == prefix,]$dataset,
            "...", fill=TRUE)
        # load the corresponding biomart dataset
        return(useMart(biomart="ensembl",
                       dataset=prefix2dataset[prefix2dataset$prefix == prefix,]$dataset))
      }
      # Otherwise return an error and stop
      else{
        stop(prefix, " is not recognised as an ensembl gene id prefix.")
      }
    }
    # If the gene id starts with "WBgene"
    else if (length(grep(pattern="^WBGene", x=sample_gene))) {
      cat("Looks like you need ensembl dataset: celegans_gene_ensembl ...", fill=TRUE)
      # load the corresponding biomart dataset
      return(useMart(biomart="ensembl", dataset="celegans_gene_ensembl"))
    }
    # If the gene id starts with "FBgn"
    else if (length(grep(pattern="^FBgn", x=sample_gene))) {
      cat("Looks like you need ensembl dataset: dmelanogaster_gene_ensembl ...", fill=TRUE)
      # load the corresponding biomart dataset
      return(useMart(biomart="ensembl", dataset="dmelanogaster_gene_ensembl"))
    }
    # If the gene id starts with "Y"
    else if (length(grep(pattern="^Y", x=sample_gene))) {
      cat("Looks like you need ensembl dataset: scerevisiae_gene_ensembl ...", fill=TRUE)
      # load the corresponding biomart dataset
      return(useMart(biomart="ensembl", dataset="scerevisiae_gene_ensembl"))
    }
    # If the gene id does not match any known ensembl gene id prefix, return an error and stop
    else{
      stop(sample_gene, " is not recognised as an ensembl gene id. Cannot guess the species.")
    }
  }
}