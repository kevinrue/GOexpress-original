cluster_GO = function(
    go_id, result, expr_data, phenodata, f=result$factor, 
    method_dist="euclidean", method_hclust="average", cex=0.8,
    main=paste(go_id, result$GO[result$GO$go_id == go_id, "name_1006"]),
    xlab="Distance", ...){
    # Fetch the list of genes associated with the go_id
    gene_ids = list_genes(go_id=go_id, result=result, data.only=TRUE)
    # Fetch and transform the expression data for those genes
    genes_expr = t(expr_data[gene_ids,])
    # Hierarchical clustering using dist and hclust
    # (The clearest method to read the labels and control their size)
    di <- dist(genes_expr, method=method_dist, ...)
    cl <- hclust(di, method=method_hclust, ...)
    # Rows are samples, label them according to the user's chosen factor
    sample_labels = Biobase::pData(phenodata)[,f]
    plot(cl, hang=-1, label=sample_labels, cex=cex, main=main, xlab=xlab, ...)
}

expression_plot = function(
    gene_id, result, expr_data, phenodata, x_var, 
    f=result$factor, ylab = "log2(cpm)", col.palette="Accent",
    col = RColorBrewer::brewer.pal(n=length(levels(Biobase::pData(
        phenodata)[,f])), name=col.palette), level=0.95, title=NULL,
    title.size=2, axis.title.size=20, axis.text.size=15,
    legend.text.size=15, legend.title.size=20, legend.key.size=30){
    # if the gene identifier is absent from the dataset
    if (!gene_id %in% rownames(expr_data)){
        # suggest close matches if any
        matches = agrep(pattern=gene_id, x=rownames(expr_data),
                        max.distance = 1, fixed=TRUE, value=TRUE)
        if (length(matches) > 0){
            cat(gene_id, "not found in dataset. Did you mean:", fill=TRUE)
            return(matches)
        }
        else{
            cat(gene_id, "not found in dataset. No close match either.")
            return()
        }
    }
    # if the result variable provided does not contain the essential slot
    if (! "genes" %in% names(result)){
        # return an error and stop
        stop("\"result=\" argument does not look like a GO_analyse output.")
    }
    # If the X variable requested does not exist in the sample annotations
    if (! x_var %in% colnames(Biobase::pData(phenodata))){
        # Return an error and stop
        stop("\"x_var=\" argument is not a valid factor in pData(phenodata).")
    }
    # Build the title message from the combination of gene_id and gene_symbol
    title = paste(gene_id, " = ", result$genes[gene_id,]$external_gene_id)
    # Assemble a data frame containing the necessary information for ggplot
    df = data.frame(Expression=expr_data[gene_id,],
                    Factor=Biobase::pData(phenodata)[,f],
                    X=Biobase::pData(phenodata)[,x_var])
    # Generate the plot
    gg = ggplot2::ggplot(df) +
        ggplot2::geom_smooth(aes(x=X, y=Expression, group = Factor,
                                 color = Factor, fill=Factor), level=level) +
        ggplot2::labs(title=title, x=x_var, y=ylab) +
        ggplot2::theme(
            plot.title = ggplot2::element_text(
                            size = ggplot2::rel(title.size)),
            axis.title=ggplot2::element_text(size=axis.title.size),
            axis.text=ggplot2::element_text(size=axis.text.size),
            legend.text=ggplot2::element_text(size=legend.text.size),
            legend.title=ggplot2::element_text(size=legend.title.size),
            legend.key.size=unit(legend.key.size, "points")) +
        ggplot2::scale_colour_manual(values=col, name=f) + 
        ggplot2::scale_fill_manual(values=col, name=f)
    # Return the plot
    return(gg)
}


expression_plot_symbol = function(
    gene_symbol, result, expr_data, phenodata, x_var, f=result$factor,
    index=0, ylab="log2cpm", col.palette="Accent",
    col = RColorBrewer::brewer.pal(n=length(levels(Biobase::pData(
        phenodata)[,f])), name=col.palette), level=0.95, titles=c(),
    title.size=2, axis.title.size=20, axis.text.size=15,
    legend.text.size=20, legend.title.size=20, legend.key.size=30){
    # if the result provided does not look like it should
    if (! "genes" %in% names(result)){
        stop("\"result=\" argument does not look like a GO_analyse output.")
    }
    # If the X variable requested does not exist in the sample annotations
    if (! x_var %in% colnames(Biobase::pData(phenodata))){
        stop("\"x_var=\" argument is not a valid factor in pData(phenodata).")
    }
    # the GO_analyse result provided contains the annotation of each gene
    # identifier
    # present in the dataset to a gene name, if any
    cat("Fetching ensembl identifier(s) annotated to",    gene_symbol, "...",
        fill=TRUE)
    mapping = data.frame(gene_id=rownames(result$genes), 
                         external_gene_id=result$genes$external_gene_id,
                         stringsAsFactors=FALSE)
    # if the gene name is absent from the mapping table
    if(!gene_symbol %in% mapping$external_gene_id){
        # suggest close matches if any
        matches = agrep(pattern=gene_symbol, x=mapping$external_gene_id,
                        fixed=TRUE, value=TRUE)
        # if we do have one or more close matches to the symbol
        if (length(matches) > 0){
            # list them to the user for help and stop the function
            cat(gene_symbol, "not found in dataset. Did you mean:", fill=TRUE)
            return(matches)
        }
        # if we don't have close matches in the dataset, tell the user and stop
        # the function
        else{
            stop(paste(gene_symbol, "not found in dataset. No close match either."))
        }
    }
    # At this stage we know the gene symbol has at least one corresponding gene
    # identifier in the BioMart database, fetch all identifier(s) corresponding
    # to that gene symbol
    gene_ids = mapping[mapping$external_gene_id == gene_symbol, "gene_id"]
    # However, we still don't know how many of those identifiers are present in
    # the expression dataset. Remove the ensembl identifiers absent from our
    # dataset as we cannot plot them
    gene_ids_present = gene_ids[gene_ids %in% rownames(expr_data)]
    # If none of the ensembl identifiers are present in the dataset
    if (length(gene_ids_present) == 0){
        cat("Gene identifiers were found for", gene_symbol, "\n",
            "but none of them were found in the expression dataset.\n",
            "Gene identifiers were:")
        return(gene_ids)
    }
    # At this stage we are finally sure that at least one of the ensembl
    # identifiers corresponding to the gene symbol are also present in the
    # expression dataset. This is the best moment to generate as many titles as
    # there are ensembl identifier(s) annotated to the given gene symbol
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
            stop("The number of titles (", length(titles),
                 ") does not match the number of plots (",
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
            cat("Use argument \"index=1\" to plot the first gene id alone,",
                "and so on.", fill=TRUE)
            # Prepare a grid to plot multiple graphs while optimising the
            # number of columns and rows
            columns = ceiling(sqrt(length(gene_ids_present)))
            # Store all the plots in a list
            plots <- list()
            for (i in seq(1,length(gene_ids_present))){
                cat("Plotting", gene_ids_present[i], fill=TRUE)
                plots[[i]] = expression_plot(
                    gene_id=gene_ids_present[i],
                    expr_data=expr_data, phenodata=phenodata,
                    x_var=x_var, result=result, f=f, ylab = "log2(cpm)",
                    col.palette="Accent", col=col, level=level,
                    title=titles[i], title.size=title.size,
                    axis.title.size=axis.title.size,
                    axis.text.size=axis.text.size,
                    legend.text.size=legend.text.size,
                    legend.title.size=legend.title.size,
                    legend.key.size=legend.key.size)
            }
            # Plot all the graphs in the optimised lattice, using the
            # ensembl-based plotting function
            multiplot(plotlist = plots, cols = columns)
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
                cat("Plotting", gene_ids_present[index], fill=TRUE)
                expression_plot(
                    gene_id=gene_ids_present[index],
                    expr_data=expr_data, phenodata=phenodata, x_var=x_var,
                    result=result, f=f, ylab = "log2(cpm)",
                    col.palette="Accent", col=col, level=level,
                    title=titles[index], title.size=title.size,
                    axis.title.size=axis.title.size,
                    axis.text.size=axis.text.size,
                    legend.text.size=legend.text.size,
                    legend.title.size=legend.title.size,
                    legend.key.size=legend.key.size)
            }
        }
    }
    # If there is a unique gene id associated to the gene symbol
    else{
        cat("Unique gene id found for", gene_symbol, fill=TRUE)
        cat("Plotting", gene_ids_present, fill=TRUE)
        expression_plot(
            gene_id=gene_ids_present,
            expr_data=expr_data, phenodata=phenodata, x_var=x_var, 
            result=result, f=f, ylab = ylab, col.palette="Accent", col = col,
            level=level, title=titles, title.size=title.size,
            axis.title.size=axis.title.size,
            axis.text.size=axis.text.size,
            legend.text.size=legend.text.size,
            legend.title.size=legend.title.size,
            legend.key.size=legend.key.size)
    }
}

heatmap_GO = function(
    go_id, result, expr_data, phenodata, gene_names=TRUE,
    f=result$factor, scale="none", cexCol=1.2, cexRow=0.5, 
    trace="none", expr.col=gplots::bluered(75), 
    row.col.palette="Accent",
    row.col=brewer.pal(n=length(unique(Biobase::pData(
        phenodata)[,f])), name=row.col.palette),
    main=paste(go_id, result$GO[result$GO$go_id == go_id,
                                "name_1006"]),
    ...){
    # Fetch the list of genes associated with the go_id
    gene_ids = list_genes(go_id=go_id, result=result, data.only=TRUE)
    # Fetch and format the expression data for those genes
    genes_expr = t(expr_data[gene_ids,])
    # Rows are samples, label them according to the user's choson factor
    sample_labels = pData(phenodata)[,f]
    # Columns are genes, label them by identifier or name
    if (gene_names){
        gene_labels = result$genes[gene_ids,]$external_gene_id
    }
    else{
        gene_labels = gene_ids
    }
    # A vector detailing the color of each sample must be prepared
    samples.col = row.col[as.factor(pData(phenodata)[,f])]
    # Plot the heatmap of the data
    heatmap.2(genes_expr, labRow=sample_labels, labCol=gene_labels,
              scale=scale, cexCol=cexCol, cexRow=cexRow, main=main,
              trace=trace, RowSideColors=samples.col, col=expr.col, ...)
}

hist_scores = function(
    result,
    main=paste("Distribution of average scores in",
               deparse(substitute(result))), xlab="Average score", ...){
    hist(result$GO$ave_score, main=main, xlab=xlab, ...)
}

list_genes = function(go_id, result, data.only=TRUE){
    # If the go_id requested is not present in the dataset
    if (!go_id %in% result$mapping$go_id){
        # return an error and stop
        stop(go_id, "is not a valid go_id in the dataset.")
    }
    # List of gene_ids annotated to the GO term
    gene_ids = result$mapping[result$mapping$go_id == go_id,"gene_id"]
    # BioMart can have more genes associated with the GO term
    # than the number actually present in the dataset
    if (data.only){
        # Filter out those not in the dataset
        gene_ids = gene_ids[gene_ids %in% rownames(result$genes)]
    }
    # return the list of gene_ids associated with it
    return(gene_ids)
}

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
    # source:
    # http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)
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
        pushViewport(viewport(layout = grid.layout(nrow(layout),
                                                   ncol(layout))))
        
        # Make each plot, in the correct location
        for (i in 1:numPlots) {
            # Get the i,j matrix positions of the regions that contain this
            # subplot
            matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
            
            print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                            layout.pos.col = matchidx$col))
        }
    }
}

overlap_GO = function(go_ids, result, filename, mar=rep(0.1, 4), ...){
    # Check that all GO terms are present in the result variable
    if (!all(go_ids %in% result$GO$go_id)){
        stop("Some go_id(s) are absent from the result variable.")
    }
    # Check that 2-5 GO terms are given, otherwise venn.diagram would crash
    if (length(go_ids) > 5){
        stop("venn.diagram() supports at most 5 groups: Too many given.")
    }
    if (length(go_ids) < 2){
        stop("venn.diagram() requires at least 2 groups: Too few given.")
    }
    # Creates a list of 2-5 gene lists to compare
    gene_sets = list()
    for (index in 1:length(go_ids)){
        gene_sets[[index]] = list_genes(go_id=go_ids[[index]], result=result)
    }
    # Print the venn diagram to the filename
    VennDiagram::venn.diagram(x=gene_sets,
                              filename=filename,
                              category.names = go_ids,
                              mar=mar, ...)
}

plot_design = function(
    go_id, result, expr_data, phenodata,
    factors=colnames(Biobase::pData(phenodata)), main="", ...){
    # if the user changed the default value
    # check that all given factors exist in colnames(phenodata)
    if (any(factors != colnames(Biobase::pData(phenodata)))){
        for(f in factors){
            if (!f %in% colnames(Biobase::pData(phenodata))){
                # Otherwise, stop and return which of them is not a valid
                # factor name
                stop(f, " is not a valid factor name in colnames(phenodata).")
            }
        } 
    }
    # Fetch the list of genes associated with the go_id
    # If the user gave the output of a GO_analyse command as result=
    # that list contains the mapping between ensembl genes and GO_id
    gene_ids_present = list_genes(go_id=go_id, result=result,
                                             data.only=TRUE)
    GO_name = result$GO[result$GO$go_id == go_id, "name_1006"]
    # Prepare a temporary data frame plot.design-friendly
    df = data.frame(t(expr_data[gene_ids_present,]),
                    Biobase::pData(phenodata[,factors]))
    # If no custom title was given
    if (main == ""){
        # Generate a smart one (careful: the same title will be used for all
        # genes in the GO term)
        # Smart title is the name of the GO term
        main=paste(go_id, GO_name)
    }
    # Perform a plot.design of all the genes in the data frame (= in the GO
    # term and in the dataset)
    plot.design(df, main=main, ...)
}

quantiles_scores = function(result, probs=c(0.9, 0.95, 0.99, 0.999, 0.9999),
                            quartiles=FALSE){
    # If user changes to "quartiles=TRUE" then the following will be true
    if (quartiles){
        quantile(x=result$GO$ave_score)
    }
    # Default is to give range of top 10%, 5%, 1%, 0.1% and 0.01%
    else {
        quantile(x=result$GO$ave_score, probs=probs)
    }
}

rerank = function(result, rank.by="rank"){
    # Reorder the GO and gene tables accordin to the user's choice
    if (rank.by == "rank"){
        result$GO = result$GO[order(result$GO$ave_rank),]
        result$genes = result$genes[order(result$genes$Rank),]
    }
    else if (rank.by == "score") {
        result$GO = result$GO[order(result$GO$ave_score, decreasing=TRUE),]
        result$genes = result$genes[order(
            result$genes$Score, decreasing=TRUE),]
    }
    else{
        stop("Invalid ranking method: ", rank.by)
    }
    return(result)
}

subset_scores = function(result, ...){
    # Save the list of filter and value for easier referencing
    filters = list(...)
    # prepares a table where the filtering results will be saved
    filtered = data.frame(row.names=result$GO$go_id)
    # For each filter
    for (filter in names(list(...))){
        # Save the filter status of each row for this filter
        ## Filter on the total count of genes associated with the GO term
        if (filter %in% c("total_count", "total")){
            #cat(filter, "equal or more than", filters[[filter]], fill=TRUE)
            filtered[,filter] = result$GO[,"total_count"] >= filters[filter]
        }
        ## Filter on the count of genes in the dataset associated with the GO
        ## term
        else if (filter %in% c("data_count", "data")){
            #cat(filter, "equal or more than", filters[[filter]], fill=TRUE)
            filtered[,filter] = result$GO[,"data_count"] >= filters[filter]
        }
        ## Filters on the average rank of the genes associated to the GO term
        else if (filter %in% c("ave_rank")){
            #cat(filter, "equal or lower than", filters[[filter]], fill=TRUE)
            filtered[,filter] = result$GO[,filter] <= filters[filter]
        }
        ## Filters on the average score of the genes associated to the GO term
        else if (filter %in% c("ave_score")){
            #cat(filter, "equal or lower than", filters[[filter]], fill=TRUE)
            filtered[,filter] = result$GO[,filter] >= filters[filter]
        }
        ## Filters on the namespace of the GO term
        else if (filter %in% c("namespace_1003", "namespace")){
            #cat(filter, "equal to", filters[[filter]], fill=TRUE)
            # GO namespace filtering should offer shortcuts
            if (filters[filter] %in% c("biological_process", "BP")){
                filtered[,filter] = result$GO$namespace_1003 ==
                    "biological_process"
            }
            else if (filters[filter] %in% c("molecular_function", "MF")){
                filtered[,filter] = result$GO$namespace_1003 ==
                    "molecular_function"
            }
            else if (filters[filter] %in% c("cellular_component", "CC")){
                filtered[,filter] = result$GO$namespace_1003 ==
                    "cellular_component"
            }
            else{
                stop("Valid values for filter ", filter, " are ",
                     "\"biological_process\", \"BP\",",
                     "\"molecular_function\", \"MF\",",
                     "\"cellular_component\", and\"CC\".")
            }
        }
        ## other filter names cause an error
        else{
            stop(filter, " is not a valid filter.")
        }
    }
    # Only the rows passing all filters will be kept
    # (all test summarised in one yes if all test are passed)
    filtered$merge = apply(X=filtered, MARGIN=1, FUN=all)
    # Filter the different slots of results to remain consistent between them
    # and save memory
    ## Subset the score table according to the above filters
    result$GO = result$GO[filtered$merge,]
    ## Subset the gene/GO mapping to keep only the GO terms left in the score
    ## table
    result$mapping = result$mapping[result$mapping$go_id %in% result$GO$go_id,]
    ## Subset the anova table to keep only the genes annotated to the genes
    ## left in the mapping table
    result$genes = result$genes[rownames(result$genes) %in%
                                    result$mapping$gene_id,]
    return(result)
}

table_genes = function(go_id, result, data.only=FALSE){
    # If the go_id requested is not present in the dataset
    if (!go_id %in% result$mapping$go_id){
        # return an error and stop
        stop(go_id, "is not a valid go_id in the dataset.")
    }
    # Otherwise fetch the list of gene_ids associated with it
    gene_ids = list_genes(go_id=go_id, result=result, data.only=data.only)
    # Then fetch the results of those genes
    res_table = result$genes[gene_ids,]
    # If gene absent from dataset were also requested (data.only=FALSE)
    # Those do not have analysis results and therefore return NA
    # Leave NAs for the results, but put their name again as the row name
    rownames(res_table) = gene_ids
    # Return the information for all those genes
    return(res_table)
}
