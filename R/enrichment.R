#' Perform Enrichment Analysis for Key Events
#'
#' Conduct enrichment analysis for key events (KEs) using the provided data and parameters.
#'
#' @param GList A named list of data frames. Each element of the list represents
#' a specific experimental condition. The dataframe must have a column with gene ids named 'Feature', and optional numerical
#' values associated to the genes (e.g. t-test statistics, fold-changes, bmd values, etc)
#' @param list_gene_sets One of the list of gene sets for enrichment analysis provided by the package
#' @param only_significant Logical indicating whether to include only significant results.
#' @param pval_th The p-value threshold for significance.
#' @param adj.method The method for p-value adjustment.
#' @param merge_by The variable for merging the results.
#' @param numerical_properties vector of strings, corresponding to colnames in the dataframes of the GList list that indicates numerical properties of the Features
#' @param background A vector of background genes to compare against.
#' @param verbose, boolean if TRUE a txtProgressBar is printed
#' @param aggregation_function Function to be used to aggregate the numerical values of the KEs/AOPs. Possible values are median, mean, min and quantile_05. When using quantile_05 the 5% quantile of the distribution of the values is computed.
#' @return A list of enriched key events (KEs) and associated results.
#' @export
enrich_KEs_AOPs = function(GList,
                           list_gene_sets,
                           only_significant = FALSE,
                           pval_th = 0.05,
                           adj.method = "fdr",
                           merge_by = "Ke",
                           numerical_properties="Score",
                           background = unique(unlist(list_gene_sets)),
                           verbose = TRUE,
                           aggregation_function = "median"){

  # The background are all the genes that are mapped to KEs/AOPs
  # background = unique(unlist(list_gene_sets))

  # Gene in the annotated gene sets (KE/AOP) are filtered for those not present in the background
  for(i in 1:length(list_gene_sets)){
    gi = list_gene_sets[[i]]
    gi = gi[gi %in% background]
    list_gene_sets[[i]] = gi
  }

  # The assumption is that all the genes in input are in the background. However, if this does not happen, they are removed
  for(i in 1:length(GList)){
    gi = GList[[i]]
    gi = gi[gi[,"Feature"] %in% background,]
    GList[[i]] = gi
  }

  ke_enrichment_results = loop_enrichment_v2(GList = GList,
                                             list_gene_sets = list_gene_sets,
                                             background = background,
                                             only_significant = only_significant,
                                             pval_th = pval_th,
                                             adj.method = adj.method,
                                             merge_by = merge_by,numerical_properties,
                                             verbose=verbose,
                                             aggregation_function=aggregation_function)

  return(ke_enrichment_results)
}

# Perform Enrichment Analysis
#
# This function performs enrichment analysis using Fisher's exact test for gene sets based on provided statistical data and gene sets.
#
# @param genes A vector of genesets
# @param reference A vector of reference genes (the background).
# @param genesets A list of gene sets.
# @param adj The method used for adjusting p-values.
# @param verbose Logical, whether to display progress messages.
#
# @return A data frame with enriched gene sets and associated statistical information.
#
# @export
# geneset_enrichment = function(genes, reference, genesets, adj = "fdr", verbose = FALSE)
# {
#
#   tab = lapply(1:length(genesets), function(i) {
#     if (verbose == TRUE) {
#       cat("processing term", i, names(genesets)[i], "\n")
#     }
#     reference = reference[!reference %in% genes]
#     RinSet = sum(reference %in% genesets[[i]])
#     RninSet = length(reference) - RinSet
#
#     genes_in_gene_set = genes[genes %in% genesets[[i]]]
#
#
#     GinSet = sum(genes %in% genesets[[i]])
#     GninSet = length(genes) - GinSet
#     fmat = matrix(c(GinSet, RinSet, GninSet, RninSet), nrow = 2, ncol = 2, byrow = F)
#     colnames(fmat) = c("inSet", "ninSet")
#     rownames(fmat) = c("genes", "reference")
#     fish = fisher.test(fmat, alternative = "greater")
#     pval = fish$p.value
#     inSet = RinSet + GinSet
#     res = c(GinSet, inSet, pval, paste(genes_in_gene_set,collapse = ";"))
#     res
#   })
#   rtab = do.call("rbind", tab)
#   rtab = data.frame(as.vector(names(genesets)), rtab)
#   rtab = rtab[order(rtab[, 4]), ]
#   colnames(rtab) = c("TermID", "genes", "all", "pval", "Genes")
#   rtab$pval = as.numeric(as.vector(rtab$pval))
#   rtab$genes = as.numeric(as.vector(rtab$genes))
#   rtab$all = as.numeric(as.vector(rtab$all))
#
#   padj = p.adjust(rtab[, 4], method = adj)
#   tab.out = data.frame(rtab, padj)
#
#   return(tab.out)
# }

#' Gene Set Enrichment Analysis
#'
#' This function performs a gene set enrichment analysis, comparing a list of genes to a reference set across multiple gene sets.
#'
#' @param background A vector of background genes to compare against.
#' @param genesets A named list of gene sets, where each element is a vector of genes belonging to that set.
#' @param adj.method Method for adjusting p-values. Default is "fdr" (False Discovery Rate).
#' @param verbose Logical, indicating whether to print progress messages. Default is FALSE.
#' @param numerical_properties A vector of column names from the genes data frame representing numerical properties to be averaged. Default includes "BMD", "BMDU", and "BMDL".
#' @param aggregation_function Function to be used to aggregate the numerical values of the KEs/AOPs. Possible values are median, mean, min and quantile_05. When using quantile_05 the 5% quantile of the distribution of the values is computed.

#' @export
geneset_enrichment = function(genes_of_interest,
                              background, genesets,
                              adj.method = "fdr",
                              verbose = FALSE,
                              numerical_properties = c("BMD","BMDU","BMDL"),
                              aggregation_function = "median"){

  print("my gene set")
  tab = lapply(1:length(genesets), function(i) {
    if (verbose == TRUE) {
      cat("processing term", i, names(genesets)[i], "\n")
    }

    genes_of_interest_in_gene_set = genes_of_interest[genes_of_interest[,"Feature"] %in% genesets[[i]],]

    res = NULL
    if(nrow(genes_of_interest_in_gene_set)>0){
      # print(i)
      if (!is.null(numerical_properties)) {
        if (nrow(genes_of_interest_in_gene_set) > 0) {
          avg_properties = c()
          for (pi in numerical_properties) {
            if(aggregation_function == "median"){
              value_aggregated = median(genes_of_interest_in_gene_set[!is.na(genes_of_interest_in_gene_set[,pi]),pi])

            }
            if(aggregation_function == "mean"){
              value_aggregated = min(genes_of_interest_in_gene_set[!is.na(genes_of_interest_in_gene_set[,pi]),pi])

            }
            if(aggregation_function == "min"){
              value_aggregated = min(genes_of_interest_in_gene_set[!is.na(genes_of_interest_in_gene_set[,pi]),pi])

            }
            if(aggregation_function == "quantile_05"){
              value_aggregated = quantile(genes_of_interest_in_gene_set[!is.na(genes_of_interest_in_gene_set[,pi]),pi],probs = 0.05)

            }

            avg_properties = c(avg_properties,value_aggregated)
          }
          names(avg_properties) = numerical_properties
        }else{
          avg_properties = rep(0, length(numerical_properties))
          names(avg_properties) = numerical_properties
        }
      }

      vector_of_genes_of_interest = genes_of_interest[,"Feature"]
      vector_of_genes_in_gene_set = genesets[[i]]

      n_genes_of_interest_in_gene_set = length(intersect(vector_of_genes_of_interest,vector_of_genes_in_gene_set))
      n_genes_of_interest_out_from_gene_set = length(setdiff(vector_of_genes_of_interest, vector_of_genes_in_gene_set))
      n_genes_not_interest_in_gene_set = length(setdiff(vector_of_genes_in_gene_set, vector_of_genes_of_interest))
      n_genes_not_interest_out_from_gene_set = length(setdiff(background, union(vector_of_genes_of_interest, vector_of_genes_in_gene_set)))

      fmat = matrix(c(n_genes_of_interest_in_gene_set, n_genes_of_interest_out_from_gene_set, n_genes_not_interest_in_gene_set, n_genes_not_interest_out_from_gene_set), nrow = 2, ncol = 2, byrow = T)
      rownames(fmat) = c("relevantGene","nonRelevantGene")
      colnames(fmat) = c("inSet","outSet")

      res = fisher.test(fmat, alternative = "greater")
      pval = res$p.value

      if(!is.null(numerical_properties)){
        res = c(names(genesets)[i], n_genes_of_interest_in_gene_set, length(vector_of_genes_in_gene_set), pval, paste(intersect(vector_of_genes_of_interest,vector_of_genes_in_gene_set),collapse = ";"), avg_properties)
      }else{
        res = c(names(genesets)[i],n_genes_of_interest_in_gene_set, length(vector_of_genes_in_gene_set), pval, paste(intersect(vector_of_genes_of_interest,vector_of_genes_in_gene_set),collapse = ";"))
      }
    }else{
      if(!is.null(numerical_properties)){
        vector_of_genes_in_gene_set = genesets[[i]]
        res = c(names(genesets)[i], 0, length(vector_of_genes_in_gene_set), 1,"", rep(0, length(numerical_properties)))
      }else{
        vector_of_genes_in_gene_set = genesets[[i]]
        res = c(names(genesets)[i],0, length(vector_of_genes_in_gene_set), 1, "")
      }
    }

    res
  })
  rtab = do.call("rbind", tab)
  rtab = rtab[order(rtab[, 4]), ]

  if(!is.null(numerical_properties)){
    colnames(rtab) = c("TermID", "relevantGenesInGeneSet", "GeneSetSize", "pval", "Genes",paste(numerical_properties))
  }else{
    colnames(rtab) = c("TermID", "relevantGenesInGeneSet", "GeneSetSize", "pval", "Genes")
  }

  rtab = as.data.frame(rtab)
  rtab$pval = as.numeric(as.vector(rtab$pval))
  rtab$relevantGenesInGeneSet = as.numeric(as.vector(rtab$relevantGenesInGeneSet))
  rtab$GeneSetSize = as.numeric(as.vector(rtab$GeneSetSize))

  if (!is.null(numerical_properties)) {
    for (pi in numerical_properties) {
      rtab[[pi]] = as.numeric(as.vector(rtab[[pi]]))
    }
  }

  idx = which(rtab$relevantGenesInGeneSet == 0)

  if(length(idx)>0){
    non_tested = rtab[idx,]
    padj = 1
    non_tested = cbind(non_tested, padj)
    rtab = rtab[-idx,]
  }else{
    non_tested = NULL
  }

  padj = p.adjust(rtab[, "pval"], method = adj.method)
  tab.out = data.frame(rtab, padj)

  tab.out = rbind(tab.out, non_tested)

  return(tab.out)
}
#' Enrichment Analysis for Multiple Experiment Genes
#'
#' This is a wrapper of the function geneset_enrichment to performs enrichment analysis for multiple experiments.
#'
#' @param GList A named list of data frames. Each element of the list represents
#' a specific experimental condition. The dataframe must have a column with gene ids, and optional numerical
#' values associated to the genes (e.g. t-test statistics, fold-changes, bmd values, etc)
#' @param list_gene_sets A list of gene sets for enrichment analysis.
#' @param background Background gene set for enrichment analysis.
# #' @param aop_ke_table_hure AOP-KE mapping table.
#' @param only_significant Logical, whether to include only significant results.
#' @param pval_th P-value threshold for significance.
#' @param adj.method Adjustment method for p-values.
#' @param merge_by Column used to merge results with the AOP-KE mapping table. Available options "Ke" or "Aop"
#' @param numerical_properties vector of strings, corresponding to colnames in the dataframes of the GList list that indicates numerical properties of the Features
#' @param verbose, boolean if TRUE a txtProgressBar is printed
#' @return A data frame containing the enriched results.
#' @param aggregation_function Function to be used to aggregate the numerical values of the KEs/AOPs. Possible values are median, mean, min and quantile_05. When using quantile_05 the 5% quantile of the distribution of the values is computed.

#' @export
loop_enrichment_v2 = function(GList,
                              list_gene_sets, background,
                              only_significant = TRUE,
                              pval_th = 0.05,
                              adj.method = "fdr",
                              merge_by = "Ke",
                              numerical_properties = c("Score"),
                              verbose = TRUE,
                              aggregation_function = "median"){
  lenrichment_result = c()

  if(verbose){
    pb = utils::txtProgressBar(min = 0, max = length(GList), style = 3)
  }

  for (i in 1:length(GList)) {
    genes_of_interest = GList[[i]]
    genesets = list_gene_sets

    enrichment_table = geneset_enrichment(genes_of_interest = genes_of_interest,
                                     background = background,
                                     genesets = list_gene_sets,
                                     adj.method = adj.method,
                                     verbose = F,
                                     numerical_properties,
                                     aggregation_function=aggregation_function) # adj can be changed to other methods too. Check documentation.

    if (only_significant) {
      enrichment_table = enrichment_table[enrichment_table$padj < pval_th,]
      if (nrow(enrichment_table) == 0 ) {
        next()
      }
    }

    ET = c()
    for (index in 1:nrow(enrichment_table)) {
      elem = unlist(strsplit(x = as.character(enrichment_table[index,"TermID"]),split = ";"))
      for (g in elem) {
        ET = rbind(ET, c(g, as.matrix(enrichment_table[index,2:ncol(enrichment_table)])))
      }
    }

    colnames(ET) = colnames(enrichment_table)
    ET = cbind(names(GList)[i], ET)
    colnames(ET)[1] = "Experiment"
    ET = as.data.frame(ET)

    ET$pval = as.numeric(as.vector(ET$pval))
    ET$padj = as.numeric(as.vector(ET$padj))
    ET$relevantGenesInGeneSet = as.numeric(as.vector(ET$relevantGenesInGeneSet))
    ET$GeneSetSize = as.numeric(as.vector(ET$GeneSetSize))

    # Adding the annotation to aop enrichemnt results.
    enrichment_table = merge(ET, aop_ke_table_hure, by.x = which(colnames(ET) == "TermID"), by.y = merge_by)

    lenrichment_result = rbind(lenrichment_result, enrichment_table)
    if(verbose){
      utils::setTxtProgressBar(pb,i)
    }
  }
  if(verbose){
    close(pb)
  }
  return(lenrichment_result)
}

#' Build AOP Enrichment Results for AOP Fingerprints
#'
#' This function constructs AOP enrichment results for AOP fingerprints based on specified criteria.
#'
#' @param aop_enrichment_results Enrichment results for AOPs.
#' @param ke_enrichment_results Enrichment results for KEs.
#' @param time_var Variable representing time points.
#' @param min_aop_length Minimum length of AOPs. Default is 6.
#' @param percentage_enriched_ke Percentage of enriched KEs for an AOP to be considered significant. Default is 33%.
#'
#' @return A list containing two data frames: detailed_results_only_enriched and detailed_results_all_ke_in_aop.
#'
#' @export
build_aop_for_aop_fingeprints = function(aop_enrichment_results,
                                         ke_enrichment_results,
                                         min_aop_length = 6,
                                         percentage_enriched_ke = 0.33
                                         # aop_ke_table_hure,
                                         # Annotate_AOPs
                                         ){

  detailed_results_only_enriched = c()
  detailed_results_all_ke_in_aop = c()
  unique_experiments = unique(ke_enrichment_results[,"Experiment"])

  for (exp_i in unique_experiments) {
      filtered_ke_results = ke_enrichment_results[ke_enrichment_results$Experiment == exp_i,] # All KEs enriched by a specific experiments
      individual_ke = unique(filtered_ke_results$TermID)
      filtered_aop_results = aop_enrichment_results[aop_enrichment_results$Experiment == exp_i,] # All AOP enriched by a specific experiment

      if (nrow(filtered_aop_results) == 0 || nrow(filtered_ke_results) == 0) {
        next
      }

      x = dplyr::distinct(aop_ke_table_hure[,c(1,2)]) # unique pairs of AOP-KEs
      aop_length = as.data.frame(table(x[,1])) #Number of KEs in each AOP

      keep_aop = c()
      filtered_aop_results2 = data.frame(filtered_aop_results, total_ke = NA, n.enriched_ke = NA, enriched_ke = NA, proportion = NA, stringsAsFactors = FALSE)

      for (i in 1:nrow(filtered_aop_results2)) {
        aop = as.character(filtered_aop_results2$TermID[i])
        appr_kes = unique(aop_ke_table_hure$Ke[which(aop_ke_table_hure$Aop == aop)]) # all KEs of said AOP
        enr_kes = intersect(appr_kes, individual_ke) # enriched KEs of a said AOP

        filtered_aop_results2$total_ke[i] <- length(appr_kes) # total nr of KE in a said AOP
        filtered_aop_results2$n.enriched_ke[i] <- length(enr_kes) # total nr of enriched KE in a said AOP
        filtered_aop_results2$enriched_ke[i] <- paste(enr_kes, collapse = "; ") # names of enriched KE in a said AOP

        proportion_enriched = length(intersect(appr_kes, individual_ke))/length(appr_kes) # proportion of enriched KE of a said AOP
        filtered_aop_results2$proportion[i] <- proportion_enriched

        # AOP fingerprint is: AOP significantly enriched AND min. 33% of its KEs are sign. enriched (min 2 KEs when length of the AOP is < 6). Even more strict cutoff is better.
        if (length(appr_kes) >= min_aop_length) {
          if (proportion_enriched > percentage_enriched_ke) {
            keep_aop = c(keep_aop, aop)
          }
        }else if (length(appr_kes) < min_aop_length) {
          if (length(enr_kes) > 1) {
            keep_aop = c(keep_aop, aop)
          }
        }
      }

      final_result = filtered_aop_results2[which(filtered_aop_results2$TermID %in% keep_aop),]


      XX1 = merge(final_result, unique(filtered_ke_results[,c("TermID","padj","Ke_description","Genes")]), by.x = "Ke",by.y = "TermID") # Only KEs that are enriched are listed with the AOPs where they belong
      XX = merge(final_result, unique(filtered_ke_results[,c("TermID","padj","Ke_description","Genes")]), by.x = "Ke",by.y = "TermID", all.x = TRUE) # All KEs are listed even if they are not enriched
      detailed_results_only_enriched = rbind(detailed_results_only_enriched,XX1)
      detailed_results_all_ke_in_aop = rbind(detailed_results_all_ke_in_aop,XX)
  }

  detailed_results_only_enriched = as.data.frame(detailed_results_only_enriched)
  colnames(detailed_results_only_enriched)[grep(pattern = "[a-z|A-Z].x",x = colnames(detailed_results_only_enriched))]  =
    gsub(pattern = ".x",replacement = "",x = colnames(detailed_results_only_enriched)[grep(pattern = "[a-z|A-Z].x",x = colnames(detailed_results_only_enriched))] )

  detailed_results_all_ke_in_aop = as.data.frame(detailed_results_all_ke_in_aop)
  colnames(detailed_results_all_ke_in_aop)[grep(pattern = "[a-z|A-Z].x",x = colnames(detailed_results_all_ke_in_aop))]  =
    gsub(pattern = ".x",replacement = "",x = colnames(detailed_results_all_ke_in_aop)[grep(pattern = "[a-z|A-Z].x",x = colnames(detailed_results_all_ke_in_aop))] )

  # print(dim(detailed_results_all_ke_in_aop))
  detailed_results_all_ke_in_aop = merge(detailed_results_all_ke_in_aop, Annotate_AOPs, by.x = "TermID", by.y = "AOP")
  # print(dim(detailed_results_all_ke_in_aop))

  # print(dim(detailed_results_only_enriched))
  detailed_results_only_enriched = merge(detailed_results_only_enriched, Annotate_AOPs, by.x = "TermID", by.y = "AOP")
  # print(dim(detailed_results_only_enriched))

  return(list("detailed_results_only_enriched" = detailed_results_only_enriched,
              "detailed_results_all_ke_in_aop" = detailed_results_all_ke_in_aop))
}

#' Filter DataFrame by Column Values
#'
#' Filters a data frame based on specified columns and filter values.
#'
#' @param mod_stats A data frame containing model statistics.
#' @param filter_column A character vector specifying the columns to filter.
#' @param filter_by A list of character vectors containing filter values for each column.
#'
#' @return Returns a filtered data frame if filtering conditions are met; otherwise, NULL.
#' @export
filter_df_no = function(mod_stats,
                        filter_column,
                        filter_by){

  if (!is.null(filter_column)) { # if null, no filter is applied

    if (length(filter_column) != length(filter_by)) {
      print("The number of column to filter, and the list of filtering attributes have a different length")
      return(NULL)
    }

    if (!all(filter_column %in% colnames(mod_stats))) { #if any of the specified column is not present in the dataframe NULL is returned
      print("The selected column is not available")
      return(NULL)
    }

    for (i in 1:length(filter_column)) {
      fc = filter_column[i]
      idx = which(mod_stats[,fc] %in% filter_by[[i]])
      if (length(idx) == 0) {
        return(NULL)
      }
      mod_stats = mod_stats[idx,]
    }
  }

  mod_stats
}

#' Render AOP Fingerprint Bubble Plot
#'
#' This function generates a bubble plot for rendering AOP fingerprint data.
#'
#' @param enrichement_data A data frame containing the enrichment data.
#' @param group_by Variable used for grouping the data.
#' @param group_by2 Second variable used for subgrouping the data.
#' @param time_var Column name of timepoint variable
#' @param filter_column Column for filtering the data.
#' @param filter_by Value for filtering the data.
#' @param is_group_by_numeric Logical, whether the grouping variable is numeric.
#' @param threshold_proportion Proportion threshold for filtering data.
#' @param text_cex Text size for labels.
#'
#' @return A bubble plot visualizing AOP fingerprint data.
#' @export
render_aop_fingerprint_bubble_plot = function(enrichement_data,
                                              group_AOPs = "SSbD_category",
                                              group_by,
                                              group_by2,
                                              x_axis_var,
                                              y_axis_var = "TermID",
                                              filter_column,
                                              filter_by,
                                              is_group_by_numeric,
                                              threshold_proportion,
                                              text_cex = 12){
  if (!is.null(filter_column)) {
    enrichement_data = filter_df_no(mod_stats = enrichement_data, filter_column = filter_column,filter_by = filter_by)
    if (is.null(enrichement_data)) return(make_empty_plot())
  }

  if (is_group_by_numeric) {
    enrichement_data[,group_by] = factor( enrichement_data[,group_by], levels = sort(unique(as.numeric(as.vector( enrichement_data[,group_by])))))
  }

  enrichement_data[,x_axis_var] = factor(enrichement_data[,x_axis_var], levels = sort(unique(enrichement_data[,x_axis_var])))

  if (!is.null(threshold_proportion)) {
    good_idx = which(as.numeric(as.vector(enrichement_data$proportion)) > threshold_proportion)
    if (length(good_idx) > 0) {
      enrichement_data = enrichement_data[good_idx,]
    }
  }

  p = ggplot2::ggplot(enrichement_data, ggplot2::aes(x = enrichement_data[,x_axis_var], y = enrichement_data[,y_axis_var], size = proportion, color = padj)) +
    ggplot2::geom_point(alpha = 1) +
    ggplot2::ylab("") +
    ggplot2::guides(size  = guide_legend(title = "Proportion of KEs"),
                    color = guide_colorbar(title = "AdjP", size = 1)) +
    ggplot2::xlab("") +
    # ggplot2::scale_size(range = c(2,6)) +
    ggplot2::theme_minimal(base_size = text_cex) +
    ggplot2::scale_x_discrete(position = "top") #+ facet_grid(~Experiment)

  if(!is.null(group_AOPs)){
    p = p + ggplot2::facet_grid(rows = group_AOPs,scales='free_y',space="free_y")
    # p = p + ggplot2::facet_grid(rows = "SSbD_category",scales='free_y',space="free_y")
  }

  if (!is.null(group_by)) {
    if (!is.null(group_by2)) {
      f = stats::formula(paste0(group_by2,"~",group_by,sep = ""))
    }else{
      f = stats::formula(paste0("~",group_by,sep = ""))
    }
    p = p + ggplot2::facet_grid(rows = f,scales = "free_y")
  }

  return(p)
}

#' Plot KEs average BMD by tissue
#'
#'
#' @param ke_enrichment_results A data frame containing the enrichment data.
#' @param relevant_tissues list of tissues to add in the plot.
#' @param filter_experiment a vector of strings of experiment IDs to be plotted. If NULL no filtering is performed
#' @param pheno_colnames vector of string indicating the variables from KE enrichment results to be considered for plotting,
#' @param group_by_time a string for timeid. If NULL all time points are included in the plot. Use NULL also if timepoint is not a variable inclued in the experimental setup
#'
#' @return A grouped bar plot
#' @export
group_enriched_ke_by_tissue = function(ke_enrichment_results,
                                       filter_experiment = NULL,
                                       group_by_time = NULL,
                                       pheno_colnames = c("TermID","Experiment","organ_tissue","BMD","Ke_description"),
                                       relevant_tissues = c("liver","kidney","blood", "lung",
                                                            "brain","brain/nervous tissue",
                                                            "brain/prefrontal cortex/liver")){

  if(!is.null(filter_experiment)){
    ke_enrichment_results = ke_enrichment_results[ke_enrichment_results$Experiment %in% filter_experiment,]
  }

  if(sum(relevant_tissues %in% Biological_system_annotations$organ_tissue) == length(relevant_tissues)){
    KE_annotated = merge(x = ke_enrichment_results, y = Biological_system_annotations, by.x = "TermID",by.y = "ke")

    df = KE_annotated[pheno_colnames]
    df$organ_tissue[which(is.na(df$organ_tissue))] = "uncategorized"

    df = df[df$organ_tissue %in% relevant_tissues ,]
    df$BMD = as.numeric(df$BMD)
    df = unique(df)
    df = df[order(df$BMD),]
    df$Ke_description = factor(df$Ke_description, levels = unique(df$Ke_description))

    p = ggplot(df, aes(x=Ke_description, y=BMD)) +
      geom_bar(stat="identity", fill="#f68060", alpha=.6, width=.4) +
      coord_flip() +
      xlab("") +
      theme_bw()+
      theme(text = element_text(size=15))

    if(!is.null(group_by_time)){
      p = p + facet_grid(paste0(group_by_time,"~organ_tissue",sep=""), scales = "free_y")

    }else{
      p = p + facet_grid(Experiment~organ_tissue, scales = "free_y")
    }

    return(p)
  }else{
    return(NULL)
  }

}

#' Load Ensembl Gene Annotation Data from biomaRt
#'
#' Loads gene annotation data for human, mouse, and rat from the December 2021 archive of Ensembl using the `biomaRt` package.
#' This function retrieves Ensembl gene IDs, gene symbols, Entrez IDs, and descriptions for each species.
#' The resulting data frames (`genes_human`, `genes_mouse`, `genes_rat`) are stored in the global environment.
#'
#' @return None. Side effect: creates global variables `genes_human`, `genes_mouse`, and `genes_rat`.
#' @import biomaRt
#' @export
load_biomaRt_data = function(){
  # library(biomaRt)
  # # had to use archived version
  # human = useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
  # mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
  # rat = useMart("ensembl", dataset = "rnorvegicus_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
  #
  # genes_human = getBM(attributes = c("ensembl_gene_id", "hgnc_symbol", "entrezgene_id", "description"), mart = human, useCache = FALSE)
  # genes_mouse = getBM(attributes = c("ensembl_gene_id", "mgi_symbol", "entrezgene_id",  "description"), mart = mouse, useCache = FALSE)
  # genes_rat = getBM(attributes = c("ensembl_gene_id", "rgd_symbol","entrezgene_id", "description"), mart = rat, useCache = FALSE)

  list(genes_human = genes_human, genes_mouse = genes_mouse, genes_rat = genes_rat)
}

#' Convert Gene IDs to Symbols in Enrichment Results
#'
#' Converts gene identifiers in a data frame of enrichment results to gene symbols.
#' Assumes the genes are stored in a column named `Genes` as a semicolon-separated string.
#'
#' @param enrichment_df A data frame with a `Genes` column containing semicolon-separated gene identifiers.
#' @param gene_id_type A string indicating the type of gene ID used (e.g., `"ENSEMBL"`, `"ENTREZID"`, or `"SYMBOL"`).
#' @param organism A string indicating the organism (`"human"`, `"mouse"`, or `"rat"`).
#' @param genes_human mapping file for human genes
#' @param genes_mouse mapping file for mouse genes
#' @param genes_rat mapping file for rat genes
#'
#' @return A data frame with the `Genes` column converted to gene symbols.
#' @export
convert_enrichment_genes_to_symbols = function(enrichment_df,
                                               gene_id_type,
                                               organism,
                                               genes_human, genes_mouse, genes_rat){

  for(i in 1:nrow(enrichment_df)){
    ms_genes = enrichment_df$Genes[i]
    if(length(ms_genes)>0){
      ms_genes = unlist(strsplit(x = ms_genes,split = ";"))
      ms_genes_symb = paste(convert_genes_to_symbol(genes = ms_genes,
                                                    gene_id_type = gene_id_type,
                                                    organism=organism,
                                                    genes_human=genes_human,
                                                    genes_mouse=genes_mouse,
                                                    genes_rat=genes_rat), collapse = ";")
      enrichment_df$Genes[i] = ms_genes_symb
    }
  }

  return(enrichment_df)

}

#' Convert a Vector of Gene IDs to Gene Symbols
#'
#' Converts a vector of gene identifiers (e.g., ENSEMBL or ENTREZ IDs) to gene symbols for a specified organism.
#' Uses preloaded data frames (`genes_human`, `genes_mouse`, `genes_rat`) from `load_biomaRt_data()`.
#'
#' @param genes A character vector of gene identifiers.
#' @param gene_id_type A string indicating the input gene ID type (e.g., `"ENSEMBL"`, `"ENTREZID"`, or `"SYMBOL"`).
#' @param organism A string indicating the organism (`"human"`, `"mouse"`, or `"rat"`).
#' @param genes_human mapping file for human genes
#' @param genes_mouse mapping file for mouse genes
#' @param genes_rat mapping file for rat genes
#' @return A character vector of gene symbols.
#' @export
convert_genes_to_symbol = function(genes,
                                   gene_id_type,
                                   organism,
                                   genes_human,
                                   genes_mouse,
                                   genes_rat){


  if(gene_id_type == "SYMBOL"){
    return(genes)
  }else{

    if(organism=="human"){
      if(gene_id_type=="ENSEMBL"){
        gene_map = genes_human[,c("ensembl_gene_id","hgnc_symbol")]
      }
      if(gene_id_type=="ENTREZID"){
        gene_map = genes_human[,c("entrezgene_id","hgnc_symbol")]
      }
    }

    if(organism=="rat"){
      if(gene_id_type=="ENSEMBL"){
        gene_map = genes_rat[,c("ensembl_gene_id","rgd_symbol")]
      }
      if(gene_id_type=="ENTREZID"){
        gene_map = genes_rat[,c("entrezgene_id","rgd_symbol")]
      }
    }

    if(organism=="mouse"){
      if(gene_id_type=="ENSEMBL"){
        gene_map = genes_mouse[,c("ensembl_gene_id","mgi_symbol")]
      }
      if(gene_id_type=="ENTREZID"){
        gene_map = genes_mouse[,c("entrezgene_id","mgi_symbol")]
      }
    }

    gene_symbols = c()

    for(g in genes){
      idx = which(gene_map$ensembl_gene_id %in% g)
      if(length(idx)>0){
        gene_symbols = c(gene_symbols, paste(gene_map$hgnc_symbol[idx], collapse = ";"))
      }else{
        gene_symbols = c(gene_symbols, paste(g, collapse = ";"))
      }
    }

    return(gene_symbols)
  }
}

