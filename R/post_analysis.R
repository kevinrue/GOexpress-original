hist_scores = function(result,
                       main=paste("Distribution of average F scores in",
                                  deparse(substitute(result))),
                       xlab="Average F score", ...){
  hist(result$scores$ave.F.score, main=main, xlab=xlab, ...)
}

quantiles_scores = function(result, probs=c(0.9, 0.95, 0.99, 0.999, 0.9999),
                            quartiles=c(FALSE, TRUE)){
  # If user changes to "quartiles=TRUE" then the following will be true
  if (quartiles[1]){
    quantile(x=result$scores$ave.F.score)
  }
  # Default is to give range of top 10%, 5%, 1%, 0.1% and 0.01%
  else{
    quantile(x=result$scores$ave.F.score, probs=probs)
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
    if (filter %in% c("total_count", "data_count", "sig_count",
                      "ave.F.score")){
      #cat(filter, "equal or more than", filters[[filter]], fill=TRUE)
      filtered[,filter] = result$scores[,filter] >= filters[filter]
    }
    ## Filters where the value should be equal to the given value
    else if (filter %in% c("namespace_1003")){
      #cat(filter, "equal to", filters[[filter]], fill=TRUE)
      # GO namespace filtering should offer shortcuts
      if (filters[filter] %in% c("biological_process", "BP")){
        filtered[,filter] = result$scores$namespace_1003 == "biological_process"
      }
      else if (filters[filter] %in% c("molecular_function", "MF")){
        filtered[,filter] = result$scores$namespace_1003 == "molecular_function"
      }
      else if (filters[filter] %in% c("cellular_component", "CC")){
        filtered[,filter] = result$scores$namespace_1003 == "cellular_component"
      }
      else{
        stop("Valid values for namespace_1003= are \"biological_process\", \"BP\", 
             \"molecular_function\", \"MF\", \"cellular_component\", and\"CC\".")
      }
    }
    ## other filters cause an error
    else{
      stop(filter, " is not a valid filter.")
    }
  }
  # Only the rows passing all filters will be kept
  filtered$merge = apply(X=filtered, MARGIN=1, FUN=all)
  # Filter the different slots of results to remain consistent between them
  ## Subset the score table according to the above filters
  new_scores = subset(result$score, filtered$merge)
  ## Subset the gene/GO mapping to keep only the GO terms left in the score table
  new_mapping = subset(result$mapping, result$mapping$go_id %in% new_scores$go_id)
  ## Subset the anova table to keep only the genes annotated to the genes left in the mapping table
  new_anova = subset(result$anova, rownames(result$anova) %in% new_mapping$gene_id)
  ## Return the filtered slots in a new list
  list(scores = new_scores,
       mapping = new_mapping,
       anova = new_anova,
       factor = result$factor)
}

list_genes = function(result, go_id){
  # If the go_id requested is not present in the dataset
  if (!go_id %in% result$mapping$go_id){
    # return an error and stop
    stop(go_id, "is not a valid go_id in the dataset.")
  }
  # List of gene_ids annotated to the GO term
  gene_ids = result$mapping[result$mapping$go_id == go_id,"gene_id"]
  # BioMart can have more genes associated with the GO term
  # than the number actually present in the dataset
  # Filter out those not in the dataset
  gene_ids_present = gene_ids[gene_ids %in% rownames(result$anova)]
  # Otherwise return the list of gene_ids associated with it
  return(gene_ids_present)
}

table_genes = function(result, go_id){
  # If the go_id requested is not present in the dataset
  if (!go_id %in% result$mapping$go_id){
    # return an error and stop
    stop(go_id, "is not a valid go_id in the dataset.")
  }
  # Otherwise fetch the list of gene_ids associated with it
  gene_ids = list_genes(result=result, go_id=go_id)
  # Return the information for those genes
  return(result$anova[gene_ids,])
}

heatmap_GO = function(go_id, result, expr_data, phenodata, gene_names=TRUE,
                        f=result$factor, scale="none", cexCol=1.2, cexRow=0.5, 
                        trace="none", expr.col=rev(redgreen(75)), row.col.palette="Accent",
                        row.col=brewer.pal(n=length(unique(pData(phenodata)[,f])), name=row.col.palette),
                        main=paste(go_id, result$scores[result$scores$go_id == go_id,]$name_1006),
                        ...){
  # Fetch the list of genes associated with the go_id
  gene_ids = list_genes(result, go_id)
  # Fetch and transform the expression data for those genes
  genes_expr = t(expr_data[gene_ids,])
  # Rows are samples, label them according to the user's choice
  sample_labels = Biobase::pData(phenodata)[,f]
  # Columns are genes, label them by identifier or name
  if (gene_names){
    gene_labels = result$anova[gene_ids,]$external_gene_id
  }
  else{
    gene_labels = gene_ids
  }
  # A vector detailing the color of each sample must be prepared
  samples.col = row.col[as.factor(pData(phenodata)[,f])]
  # Plot the heatmap of the data
  heatmap.2(genes_expr, labRow=sample_labels, labCol=gene_labels, scale=scale,
            cexCol=cexCol, cexRow=cexRow, main=main, trace=trace, 
            RowSideColors=samples.col, col=expr.col, ...)
}


cluster_GO = function(go_id, result, expr_data, phenodata, f=result$factor, 
                      method_dist="euclidean", method_hclust="average", cex=0.8,
                      main=paste(go_id, result$scores[result$scores$go_id == go_id,]$name_1006),
                      xlab="Distance",
                      ...){
  # Fetch the list of genes associated with the go_id
  gene_ids = list_genes(result, go_id)
  # Fetch and transform the expression data for those genes
  genes_expr = t(expr_data[gene_ids,])
  # Hierarchical clustering
  # (The clearest to read the labels and control their size)
  di <- dist(genes_expr, method=method_dist, ...) # euclidean distances between the rows
  cl <- hclust(di, method=method_hclust, ...)
  # Rows are samples, label them according to the user's choice
  sample_labels = Biobase::pData(phenodata)[,f]
  plot(cl, hang=-1, label=sample_labels, cex=cex, main=main, xlab=xlab, ...)
}

expression_plot = function(gene_id, expr_data, phenodata, x_var, result, f=result$factor, 
                           ylab = "log2(cpm)", col.palette="Accent",
                           col = brewer.pal(n=length(levels(pData(phenodata)[,f])), name=col.palette),
                           level=0.95, title=NULL, title.size=2){
  # if the gene identifier is absent from the dataset
  if (!gene_id %in% rownames(expr_data)){
    # suggest close matches if any
    matches = agrep(pattern=gene_id, x=rownames(expr_data), max.distance = 1, fixed=TRUE, value=TRUE)
    if (length(matches) > 0){
      cat(gene_id, "not found in dataset. Did you mean:", fill=TRUE)
      return(matches)
    }
    else{
      cat(gene_id, "not found in dataset. No close match either.")
      return()
    }
  }
  # if the result was provided does not look like it should
  if (class(result) != "list" || length(result) != 4){
    stop("\"result=\" argument does not look like a GO_analyse output.")
  }
  
  title = paste(gene_id, " = ", result$anova[gene_id,]$external_gene_id)
  # Assemble a data frame containing the necessary information
  df = data.frame(Expression=expr_data[gene_id,],
                  Factor=Biobase::pData(phenodata)[,f],
                  # The factor
                  X=Biobase::pData(phenodata)[,x_var]) # targets data loaded above, Animal field in it
  # Generate the plot
  gg = ggplot(df) +
    geom_smooth(aes(x=X, y=Expression, group = Factor, color = Factor, fill=Factor), level=level) +
    labs(title = title) +
    xlab(x_var) +
    ylab(ylab) +
    theme(plot.title = element_text(size = rel(title.size))) +
    scale_colour_manual(values=col) + 
    scale_fill_manual(values=col)
  # Return the plot (test)
  return(gg)
}


expression_plot_symbol = function(gene_symbol, expr_data, phenodata, x_var, result,
                                  f=result$factor, index=0, ylab="log2cpm",
                                  col.palette="Accent",
                                  col = brewer.pal(n=length(levels(pData(phenodata)[,f])), name=col.palette),
                                  level=0.95,
                                  titles=c(), title.size=2){
  # if the result was provided does not look like it should
  if (class(result) != "list" || length(result) != 4){
    stop("\"result=\" argument does not look like a GO_analyse output.")
  }
  # the GO_analyse result provided contains the annotation of each ensembl identifier
  # present in the dataset to a gene name, if any
  cat("Fetching ensembl identifier(s) annotated to",  gene_symbol, "...", fill=TRUE)
  mapping = data.frame(gene_id=rownames(result$anova), 
                       external_gene_id=result$anova$external_gene_id,
                       stringsAsFactors=F)
  # if the gene name is absent from the mapping table
  if(!gene_symbol %in% mapping$external_gene_id){
    # suggest close matches if any
    matches = agrep(pattern=gene_symbol, x=mapping$external_gene_id, fixed=TRUE, value=TRUE)
    # if we do have one or more close matches to the symbol
    if (length(matches) > 0){
      # list them to the user for help and stop the function
      cat(gene_symbol, "not found in dataset. Did you mean:", fill=TRUE)
      return(matches)
    }
    # if we don't have close matches in the dataset, tell the user and stop the function
    else{
      return(gene_symbol, "not found in dataset. No close match either.")
    }
  }
  # At this stage we know the gene symbol has at least one corresponding ensembl identifier
  # in the BioMart database, fetch all identifier(s) corresponding to that gene symbol
  gene_ids = mapping[mapping$external_gene_id == gene_symbol, "gene_id"]
  # However, we still don't know how many of those are present in the expression
  # dataset. Remove the ensembl identifiers absent from our dataset
  gene_ids_present = gene_ids[gene_ids %in% rownames(expr_data)]
  # If none of the ensembl identifiers are present in the dataset
  if (length(gene_ids_present) == 0){
    cat("Gene identifiers were found for", gene_symbol, "\n",
         "but none of them were found in the expression dataset.\n",
         "Gene identifiers were:")
    return(gene_ids)
  }
  # At this stage we are finally sure that at least one of the ensembl identifiers 
  # corresponding to the gene symbol are also present in the expression dataset
  # This is the best moment to generate as many titles as there are ensembl identifier(s)
  # annotated to the given gene symbol
  # If the user left the default vector of titles (empty)
  if (is.null(titles)){
    # Create a smart title for each plot
    for (ensembl in gene_ids_present){
      titles = c(titles, paste(gene_symbol, " = ", ensembl))
    }
  }
  # If the user changed the default vector of titles
  else{
    # if the number of titles does not match the number of plots
    if(length(titles) != length(gene_ids_present)){
      # return an error and stop
      stop("The number of titles (", length(titles), ") does not match the number of plots (",
           length(gene_ids_present), ").")
    }
  }
  # If there are strictly more than 1 gene id associated with the gene symbol
  if (length(gene_ids_present) > 1){
    # Tell the user
    cat("Multiple gene ids found for", gene_symbol, fill=TRUE)
    cat("Indices are:", fill=TRUE)
    print(gene_ids_present)
    # if the user did not change the default index value (0)
    # the function will plot all ensembl ids in a lattice
    if (index==0){
      # A first time user might not know that how to plot a single plot
      cat("Use argument \"index=1\" to plot the first gene id alone, and so on.", fill=TRUE)
      # Prepare a grid to plot multiple graphs while optimising the number of columns and rows
      columns = ceiling(sqrt(length(gene_ids_present)))
      # Store all the plots in a list
      plots <- list()
      for (i in seq(1,length(gene_ids_present))){
        cat("Plotting ", gene_ids_present[i], fill=TRUE)
        plots[[i]] = expression_plot(gene_id=gene_ids_present[i], expr_data=expr_data, phenodata=phenodata,
                                     x_var=x_var, result=result, f=f, col=col,
                                     level=level, title=titles[i], title.size=title.size)
      }
      # Plot all the graphs in the optimised lattice, using the ensembl-based plotting function
      multiplot(plotlist = plots, cols = columns)
      # Restore the original lattice
      #par = opar
    }
    # If the user gave a non-zero index
    else{
      # If the index is out of bound
      if (abs(index) > length(gene_ids_present)){
        # Return an error
        print("Index is out of bound.")
        cat("Indices are:", fill=TRUE)
      }
      # If the index is acceptable
      else{
        # Plot the corresponding graph
        cat("Plotting ", gene_ids_present[index], fill=TRUE)
        expression_plot(gene_id=gene_ids_present[index], expr_data=expr_data, phenodata=phenodata,
                        x_var=x_var, result=result, f=f, col=col,
                        level=level, title=titles[index], title.size=title.size)
      }
    }
  }
  # If there is a unique gene id associated to the gene symbol
  else{
    cat("Unique gene id found for", gene_symbol, fill=TRUE)
    cat("Plotting ", gene_ids_present, fill=TRUE)
    expression_plot(gene_id=gene_ids_present, expr_data=expr_data, phenodata=phenodata, x_var=x_var,
                    result=result, f=f, ylab = ylab,
                    col = col, level=level,
                    title=titles,
                    title.size=title.size)
  }
}

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  # from http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)), byrow=TRUE,
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}


plot_design = function(go_id, expr_data, phenodata, result,
                       factors=colnames(Biobase::pData(phenodata)),
                       main="", ...){
  # if the user changed the default value
  # check that all given factors exist in colnames(phenodata)
  if (all(factors != colnames(Biobase::pData(phenodata)))){
    for(f in factors){
      if (!f %in% colnames(Biobase::pData(phenodata))){
        # Otherwise, stop and return which of them is not a valid factor name
        stop(f, " is not a valid factor name in colnames(phenodata).")
      }
    } 
  }
  # Fetch the list of genes associated with the go_id
  # if the result was provided does not look like it should
  if (class(result) != "list" || length(result) != 4){
    stop("\"result=\" argument does not look like a GO_analyse output.")
  }
  # If the user gave the output of a GO_analyse command as result=
  # that list contains the mapping between ensembl genes and GO_id
  gene_ids_present = list_genes(result, go_id)
  GO_name = result$scores[result$scores$go_id == go_id,]$name_1006
  # Prepare a temporary data frame plot.design-friendly
  df = data.frame(t(expr_data[gene_ids_present,]), Biobase::pData(phenodata[,factors]))
  # If no custom title was given
  if (main == ""){
    # Generate a smart one (careful: the same title will be used for all genes in the GO term)
    # Smart title is the name of the GO term
    main=paste(go_id, GO_name)
  }
  # Perform a plot.design of all the genes in the data frame (= in the GO term and in the dataset)
  plot.design(df, main=main, ...)
}