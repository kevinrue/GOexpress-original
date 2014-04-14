hist_scores = function(result,
                       main=paste("Distribution of average F scores in",
                                  deparse(substitute(result))),
                       xlab="Average F score", ...){
  hist(result$scores$ave.F.score, main=main, xlab=xlab, ...)
}

quantiles_scores = function(result, quartiles=c(FALSE, TRUE)){
  # If user changes to "quartiles=TRUE" then the following will be true
  if (quartiles[1]){
    quantile(x=result$scores$ave.F.score)
  }
  # Default is to give range of top 10%, 5%, 1%, 0.1% and 0.01%
  else{
    quantile(x=result$scores$ave.F.score, c(0.9, 0.95, 0.99, 0.999, 0.9999))
  }
}

subset_scores = function(result, ...){
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
      #cat(filter, "equal or more than", filters[[filter]], fill=TRUE)
      filtered[,filter] = result$scores[,filter] >= filters[filter]
    }
  }
  #
  filtered$merge = apply(X=filtered, MARGIN=1, FUN=all)
  #
  list(scores=subset(result$score, filtered$merge),
       mapping=result$mapping,
       anova=result$anova,
       factor=result$factor)
}

list_genes = function(result, go_id){
  # If the go_id requested is not present in the dataset
  if (!go_id %in% result$mapping$go_id){
    # return an error and stop
    stop(go_id, "is not a valid go_id in the dataset.")
  }
  # Otherwise return the list of ensembl_ids associated with it
  return(result$mapping[result$mapping$go_id == go_id,]$ensembl_gene_id)
}

table_genes = function(result, go_id){
  # If the go_id requested is not present in the dataset
  if (!go_id %in% result$mapping$go_id){
    # return an error and stop
    stop(go_id, "is not a valid go_id in the dataset.")
  }
  # Otherwise fetch the list of ensembl_ids associated with it
  ensembls = result$mapping[result$mapping$go_id == go_id,]$ensembl_gene_id
  # Return the information for those genes
  return(result$anova[ensembls,])
}

heatmap_GO = function(go_id, result, expr_data, phenodata, gene_names=TRUE,
                      f=result$factor, scale="none", cexCol=1.2, cexRow=0.5, 
                      main=paste(go_id, result$scores[result$scores$go_id == go_id,]$name_1006),
                      ...){
  # Fetch the list of genes associated with the go_id
  ensembl_ids = list_genes(result, go_id)
  # Fetch and transform the expression data for those genes
  genes_expr = t(expr_data[ensembl_ids,])
  # Rows are samples, label them according to the user's choice
  sample_labels = pData(phenodata)[,f]
  # Columns are genes, label them by identifier or name
  if (gene_names){
    gene_labels = result$anova[ensembl_ids,]$external_gene_id
  }
  else{
    gene_labels = ensembl_ids
  }
  # Plot the heatmap of the data
  heatmap(genes_expr, labRow=sample_labels, labCol=gene_labels, scale=scale,
          cexCol=cexCol, cexRow=cexRow, main=main, ...)
}

cluster_GO = function(go_id, result, expr_data, phenodata, f=result$factor, 
                      method_dist="euclidean", method_hclust="average", cex=0.8,
                      main=paste(go_id, result$scores[result$scores$go_id == go_id,]$name_1006),
                      xlab="Distance",
                      ...){
  # Fetch the list of genes associated with the go_id
  ensembl_ids = list_genes(result, go_id)
  # Fetch and transform the expression data for those genes
  genes_expr = t(expr_data[ensembl_ids,])
  # Hierarchical clustering
  # (The clearest to read the labels and control their size)
  di <- dist(genes_expr, method=method_dist, ...) # euclidean distances between the rows
  cl <- hclust(di, method=method_hclust, ...)
  # Rows are samples, label them according to the user's choice
  sample_labels = pData(phenodata)[,f]
  plot(cl, hang=-1, label=sample_labels, cex=cex, main=main, xlab=xlab, ...)
}