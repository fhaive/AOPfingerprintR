#' Compute closest adverse outcomes to a given vertex
#'
#' This function takes in a vertex (or set of vertices) of interest, an igraph object representing the KE-KE network,
#' a vector of adverse outcomes (AO), and optional parameters threshold and distance, and returns the closest adverse outcomes
#' to the vertex of interest in the KEKE_net. If distance is set to TRUE, it returns both the names of the closest adverse outcomes
#' and the corresponding distances between the vertex of interest and the adverse outcomes.
#'
#' @param interesting_vertex A string representing the name of the vertex/vertices of interest in the KEKE_net.
#' @param KEKE_net An igraph object representing the knowledge exchange network.
#' @param AO A vector of adverse outcomes.
#' @param max_path_length An integer representing the maximum length of the path allowed to reach an adverse outcome (Including KE and AO).
#' @param n_AOs An integer representing the maximumn number adverse outcomes to be added
#'
#' @return the function returns a vector of the names of the closest AOs to the vertex of interest in the KEKE_net.
#'
#' @export
compute_the_closest_AOs = function(interesting_vertex, KEKE_net, AO,max_path_length=10, mode = "out", n_AOs = 20) {

  from_ke = names(igraph::V(KEKE_net)[which(igraph::V(KEKE_net)$name %in% interesting_vertex)])
  to_aos =  names(igraph::V(KEKE_net)[which(igraph::V(KEKE_net)$name %in% AO)])

  res_list = list()
  for(fk in from_ke){ # for each ke of interest
    res = igraph::shortest_paths(graph = KEKE_net, from = fk,
                                 to = to_aos,
                                 mode = "all",output = "vpath") # find the sthortest path vs all Kes

    l_path = unlist(lapply(res$vpath, length)) # compute shortest pathlength
    names(l_path) = to_aos
    l_path = l_path[l_path>0]
    l_path = l_path[l_path <= max_path_length] # The path (including the starting KE and ending AO should not be longer that a certain number of KEs)
    l_path = sort(l_path,decreasing = F)

    if(length(l_path)>n_AOs){
      l_path = l_path[1:n_AOs]
    }

    idx = which(to_aos %in% names(l_path))
    res_list[[fk]] = names(unlist(res$vpath[idx]))
  }

  closest_AOs = unique(unlist(res_list))
  return(closest_AOs)
}

#' Compute closest molecular initiating events to a given vertex/set of vertices
#'
#' This function takes in a vertex (or set of vertices) of interest, an igraph object representing the KE-KE network,
#' a vector of adverse outcomes (AO), and optional parameters threshold and distance, and returns the closest adverse outcomes
#' to the vertex of interest in the KEKE_net. If distance is set to TRUE, it returns both the names of the closest adverse outcomes
#' and the corresponding distances between the vertex of interest and the adverse outcomes.
#'
#' @param interesting_vertex A string representing the name of the vertex/vertices of interest in the KEKE_net.
#' @param KEKE_net An igraph object representing the knowledge exchange network.
#' @param MIE A character vector representing the list of molecular initiating events (MIEs).
#' @param max_path_length An integer representing the maximum length of the path allowed (Including MIE and KE)
#' @param n_MIEs An integer representing the maximumn number adverse outcomes to be added#' @param distance A logical value (TRUE or FALSE) indicating whether to return both the names of the closest adverse outcomes and the corresponding distances (default is FALSE).
#'
#' @return the function returns a vector of the names of the closest MIEs to the vertex of interest in the KEKE_net.
#'
#' @export
compute_the_closest_MIEs = function(interesting_vertex, KEKE_net, MIE, max_path_length = 5, distance=TRUE, mode = "out",n_MIEs = 20){

  from_mie =  names(igraph::V(KEKE_net)[which(igraph::V(KEKE_net)$name %in% MIE)])
  to_ke = names(igraph::V(KEKE_net)[which(igraph::V(KEKE_net)$name %in% interesting_vertex)])

  res_list = list()
  for(fk in to_ke){ # for each ke of interest

    l_path = rep(0, length(from_mie))
    names(l_path) = from_mie
    for( mi in from_mie){
      res = igraph::shortest_paths(graph = KEKE_net, from = mi,
                                   to = fk,
                                   mode = "all",output = "vpath") # find the sthortest path vs all Kes
      l_path[mi] = unlist(lapply(res$vpath, length))
    }


    l_path = l_path[l_path>0]
    l_path = l_path[l_path <= max_path_length]  # The path (including the MIE and KE should not be longer that a certain number of KEs)
    l_path = sort(l_path,decreasing = F)

    if(length(l_path)>n_MIEs){
      l_path = l_path[1:n_MIEs]
    }

    idx = which(from_mie %in% names(l_path))
    res_list[[fk]] = names(unlist(res$vpath[idx]))
  }

  closest_MIEs = unique(unlist(res_list))
  return(closest_MIEs)
}

#' Plot a Visualization Network
#'
#' This function generates a visualization network (visNetwork) based on the nodes and edges derived from enrichment data. It allows for the visual grouping of nodes based on either AOPs (Adverse Outcome Pathways) or SSbD (Safe and Sustainable by Design) impact categories, and highlights key events (KEs) based on their significance.
#'
#' @param nodes A data frame containing nodes representing Key Events (KEs). Each node should include information such as KE identifiers, KE descriptions, p-values, associated genes, and any numerical variables relevant to the analysis.
#' @param edges A data frame containing edges representing the connections between the KEs in the network. Each edge should include information such as the source and target KE IDs.
#' @param numerical_variables A character vector specifying the names of the columns in `nodes` that contain numerical data to be displayed in the node tooltips.
#' @param group_by A character string indicating how to group the nodes in the network. Accepts either "aop" (default) to group by Adverse Outcome Pathways, or other specified categories such as "ssbd" to group by SSbD impact categories.Other grouping categories can be the colnames of the nodes parameter
#' @return A visNetwork object that can be rendered to create an interactive network visualization. The network will color nodes based on their significance and group them according to the specified category.
#' @export
# plot_visNetwork = function(nodes,
#                            edges,
#                            numerical_variables,
#                            group_by = "aop"){
#
#
#   color = rep("gray",nrow(nodes))
#   # color[!is.na(nodes[,numerical_variables[1]])] = "red"
#   color[which(nodes$padj<0.05)] = "red"
#   nodes$color = color
#
#   nd <- nodes$id
#   ke_type <- c()
#   for (n in 1:length(nd)) {
#     ke <- nd[n]
#     ke_annotation <- aop_ke_table_hure$Ke_type[aop_ke_table_hure$Ke == ke]
#     if (!"MolecularInitiatingEvent" %in% ke_annotation & !"AdverseOutcome" %in% ke_annotation) {
#       ke_type = c(ke_type,"KeyEvent")
#     } else {
#       tab = table(ke_annotation)
#       to_rem = which(names(tab) == "KeyEvent")
#       if (length(to_rem) > 0) tab = tab[-to_rem]
#       tab = sort(tab,decreasing = T)
#       ke_type = c(ke_type,names(tab)[1])
#     }
#   }
#
#   my_shapes = c("triangle","square","star")
#   names(my_shapes) = c("MolecularInitiatingEvent","KeyEvent","AdverseOutcome")
#   nodes$shape = my_shapes[ke_type]
#
#   aop_groups = c()
#   ssbd_indicator_groups = c()
#   for (i in nodes$id) {
#     aops = aop_ke_table_hure$a.name[aop_ke_table_hure$Ke == i]
#     if (length(aops) > 1) {
#       aops = paste(aops,collapse = ", ")
#     }
#     aop_groups = c(aop_groups, aops)
#
#     aops_id = aop_ke_table_hure$Aop[aop_ke_table_hure$Ke == i]
#     ssbd_indicator_id = unique(Annotate_AOPs$SSbD_category[Annotate_AOPs$AOP %in% aops_id])
#
#     if(length(ssbd_indicator_id) > 1){
#       ssbd_indicator_id = paste(ssbd_indicator_id, collapse = ", ")
#     }
#     ssbd_indicator_groups = c(ssbd_indicator_groups, ssbd_indicator_id)
#   }
#
#   names(aop_groups) = nodes$id
#   names(ssbd_indicator_groups) = nodes$id
#
#   if(group_by == "aop"){ # either select individual AOPs
#     nodes$group = aop_groups
#   }else{ # select ssbd impact categorids
#     nodes$group = ssbd_indicator_groups
#   }
#
#   nodes$title = paste0("<p> Id:",nodes$id, "</p>",
#                        "<p> Description: ",nodes$label, "</p>")
#
#   for(nv in numerical_variables){
#        nodes$title = paste0(nodes$title, "<p>",nv,": ", nodes[[nv]],"</p>")
#        # "<p> BMDL: ",nodes$BMDL, "</p>",
#        # "<p> BMD: ",nodes$BMD, "</p>",
#        # "<p> BMDUL ",nodes$BMDU, "</p>",
#   }
#
#   nodes$title = paste0(nodes$title, "<p> padj: ",nodes$padj, "</p>",
#                        "<p> AOPs: <p> ",aop_groups, "</p>",
#                        "<p> SSbd: <p> ",ssbd_indicator_groups, "</p>",
#                        "<p> genes: ",nodes$genes, "</p>")
#
#   vn = visNetwork::visNetwork(nodes, edges) %>%
#     visNetwork::visNodes(color = "red") %>%
#     visNetwork::visLegend(position = "right",
#                           addNodes = list(list(label = "MolecularInitiatingEvent", shape = "triangle"),
#                                           list(label = "KeyEvent", shape = "square"),
#                                           list(label = "AdverseOutcome", shape = "star"),
#                                           list(label = "Enriched", color = "red"),
#                                           list(label = "Not enriched", color = "gray")
#                           ),
#                           useGroups = FALSE) %>%
#     visNetwork::visOptions(selectedBy = list(variable = "group", multiple = TRUE))
#
#
#   return(vn)
#
# }

plot_visNetwork = function(nodes,
                           edges,
                           numerical_variables,
                           group_by = "aop",
                           color_enriched="red",
                           color_not_enriched="gray"){


  color = rep(color_not_enriched,nrow(nodes))
  # color[!is.na(nodes[,numerical_variables[1]])] = "red"
  color[which(nodes$padj<0.05)] = color_enriched
  nodes$color = color

  nd <- nodes$id
  ke_type <- c()
  for (n in 1:length(nd)) {
    ke <- nd[n]
    ke_annotation <- aop_ke_table_hure$Ke_type[aop_ke_table_hure$Ke == ke]
    if (!"MolecularInitiatingEvent" %in% ke_annotation & !"AdverseOutcome" %in% ke_annotation) { #if it is never a MIE or AO for any AOPs then is marked as KE
      ke_type = c(ke_type,"KeyEvent")
    } else {
      #if the KE is marked as MIE or AO of one or more AOP, then in the final
      #network is marked MIE or AO based on the maximum frequency of being MIE or AO across the AOPs.
      # E.g. if it is MIE twice and AO only once, is marked as MIE
      tab = table(ke_annotation)
      to_rem = which(names(tab) == "KeyEvent")
      if (length(to_rem) > 0) tab = tab[-to_rem]
      tab = sort(tab,decreasing = T)
      ke_type = c(ke_type,names(tab)[1])
    }
  }

  my_shapes = c("triangle","square","star")
  names(my_shapes) = c("MolecularInitiatingEvent","KeyEvent","AdverseOutcome")
  nodes$shape = my_shapes[ke_type]

  Annotate_AOPs$Organ[is.na(Annotate_AOPs$Organ)] = "Not annotated to organ"
  Annotate_AOPs = cbind(Annotate_AOPs,
                        "SSbD_category_organ"= paste(Annotate_AOPs$SSbD_category,
                                                    Annotate_AOPs$Organ,
                                                    sep=" - "))
  aop_groups = c()
  ssbd_indicator_groups = c()
  for (i in nodes$id) {
    aops = NULL
    aops_id = NULL
    ssbd_indicator_id = NULL

    aops = aop_ke_table_hure$a.name[aop_ke_table_hure$Ke == i]
    aops <- gsub(",", "-", aops)
    if (length(aops) > 1) {
      aops = paste(aops,collapse = ", ")
    }
    aop_groups = c(aop_groups, aops)

    aops_id = aop_ke_table_hure$Aop[aop_ke_table_hure$Ke == i]
    ssbd_indicator_id = unique(Annotate_AOPs$SSbD_category_organ[Annotate_AOPs$AOP %in% aops_id])
    ssbd_indicator_id <- gsub(",", "-", ssbd_indicator_id)

    if(length(ssbd_indicator_id) > 1){
      ssbd_indicator_id = paste(ssbd_indicator_id, collapse = ", ")
    }
    ssbd_indicator_groups = c(ssbd_indicator_groups, ssbd_indicator_id)
  }

  names(aop_groups) = nodes$id
  names(ssbd_indicator_groups) = nodes$id

  if(group_by == "aop"){ # either select individual AOPs
    nodes$group = aop_groups
  }else{
    if(group_by == "ssbd"){
      # select ssbd impact categorids
      nodes$group = ssbd_indicator_groups
    }else{
      if(group_by %in% colnames(nodes)){
        nodes$group = as.vector(nodes[,group_by])[[1]]
      }else{
        print("group_by should either be 'aop','ssbd' or another variable in colnames(nodes)")
        return(NULL)
      }
    }
  }

  nodes$title = paste0("<p> Id:",nodes$id, "</p>",
                       "<p> Description: ",nodes$label, "</p>")

  for(nv in numerical_variables){
    nodes$title = paste0(nodes$title, "<p>",nv,": ", nodes[[nv]],"</p>")
    # "<p> BMDL: ",nodes$BMDL, "</p>",
    # "<p> BMD: ",nodes$BMD, "</p>",
    # "<p> BMDUL ",nodes$BMDU, "</p>",
  }

  nodes$title = paste0(nodes$title, "<p> padj: ",nodes$padj, "</p>",
                       "<p> AOPs: <p> ",aop_groups, "</p>",
                       "<p> SSbd: <p> ",ssbd_indicator_groups, "</p>",
                       "<p> genes: ",nodes$genes, "</p>")

  vn = visNetwork::visNetwork(nodes, edges) %>%
    visNetwork::visNodes(color = color_enriched) %>%
    visNetwork::visLegend(position = "right",
                          addNodes = list(list(label = "MolecularInitiatingEvent", shape = "triangle"),
                                          list(label = "KeyEvent", shape = "square"),
                                          list(label = "AdverseOutcome", shape = "star"),
                                          list(label = "Enriched", color =color_enriched),
                                          list(label = "Not enriched", color = color_not_enriched)
                          ),
                          useGroups = FALSE) %>%
    visNetwork::visOptions(selectedBy = list(variable = "group", multiple = TRUE))


  return(vn)

}


#' Create a Visualization Network
#'
#' This function generates a visualization network (visNetwork) based on enrichment data and a network structure.
#'
#' @param detailed_results A data frame containing detailed enrichment results. It should include columns corresponding to the experiment name, KE (Key Event) identifiers, KE descriptions, and other relevant numerical variables.
#' @param experiment The name of the experiment for which the network is being created. This should match the experiment name in the `detailed_results` data frame.
#' @param enlarge_ke_selection Logical, whether to expand the selection of Key Events (KEs) to include the closest Molecular Initiating Events (MIEs) and Adverse Outcomes (AOs). Default is TRUE.
#' @param ke_id The column name in `detailed_results` that corresponds to the Key Event (KE) identifier.
#' @param numerical_variables A character vector specifying the names of the columns in `detailed_results` that contain numerical variables to be included in the network nodes.
#' @param pval_variable The name of the column in `detailed_results` that contains p-values or adjusted p-values associated with the enrichment results.
#' @param gene_variable The name of the column in `detailed_results` that contains gene identifiers or gene sets associated with the enrichment results.
#' @param convert_to_gene_symbols Logical, whether to convert gene identifiers to gene symbols using organism-specific mapping. Default is FALSE.
#' @param gene_id_type A character string indicating the type of gene identifiers used (e.g., "ENSEMBL", "ENTREZID", or "SYMBOL").
#' @param organism A character string specifying the organism of interest ("human", "mouse", or "rat").
#' @param max_path_length Maximum length of the shortest path to be used to retrieve close MIEs or AOs if `enlarge_ke_selection = TRUE`. Default is 20.
#' @param n_AOs Number of AOs to add if `enlarge_ke_selection = TRUE`. Default is 10.
#' @param n_MIEs Number of MIEs to add if `enlarge_ke_selection = TRUE`. Default is 10.
#' @param mode Parameter passed to the `igraph::shortest.paths` function if `enlarge_ke_selection = TRUE`. Options are "all", "out", or "in". Default is "out".
#' @param genes_human mapping file for human genes
#' @param genes_mouse mapping file for mouse genes
#' @param genes_rat mapping file for rat genes

#' @return A list containing two elements:
#' \describe{
#'   \item{nodes}{A data frame of nodes representing the selected KEs, including information on the numerical variables, p-values, and associated genes.}
#'   \item{edges}{A data frame of edges representing the connections between the KEs in the network.}
#' }
#' @export
make_visNetwork = function(detailed_results,
                           experiment,
                           enlarge_ke_selection = TRUE,
                           ke_id,
                           numerical_variables,
                           pval_variable,
                           gene_variable,
                           convert_to_gene_symbols = F,
                           gene_id_type,
                           organism,
                           max_path_length=20,
                           n_AOs = 10,
                           n_MIEs = 10,
                           mode = "out",
                           genes_human, genes_mouse, genes_rat){

  KEKE_net = igraph::upgrade_graph(KEKE_net)
  # retrieve list of MIEs and AOs
  MIE = unique(aop_ke_table_hure$Ke[aop_ke_table_hure$Ke_type == 'MolecularInitiatingEvent'])
  AO = unique(aop_ke_table_hure$Ke[aop_ke_table_hure$Ke_type == 'AdverseOutcome'])

  node_df = detailed_results[detailed_results$Experiment == experiment,]
  node_df = unique(node_df)

  if(nrow(node_df) == 0) return(NULL)

  interesting_vertex = unique(node_df[,ke_id])

  if (enlarge_ke_selection) {
    closest_AOs = compute_the_closest_AOs(interesting_vertex, KEKE_net, AO = AO,max_path_length = max_path_length,n_AOs = n_AOs, mode = mode)
    closest_MIEs = compute_the_closest_MIEs(interesting_vertex, KEKE_net, MIE = MIE,max_path_length = max_path_length,  n_MIEs = n_MIEs, mode = mode)
    events = union(node_df[,ke_id], union(closest_AOs, closest_MIEs))
    events = events[!is.na(events)]
    events = events[events %in% aop_ke_table_hure$Ke]
  }else{
    events = interesting_vertex
    events = events[!is.na(events)]
    events = events[events %in% aop_ke_table_hure$Ke]

  }

  library(igraph)
  # select the relevant KEs from the KEKE network
  selegoG <- igraph::induced_subgraph(KEKE_net,igraph::V(KEKE_net)[which(igraph::V(KEKE_net)$name %in% events)])

  # retrieve edge list
  edgelist <- igraph::as_edgelist(selegoG,names = T)

  if(nrow(edgelist)==0){
    return(list("nodes" = NULL,"edges" = NULL))

  }

  edges = data.frame(from = edgelist[,1], to = edgelist[,2], "title" = igraph::E(selegoG)$type.r, arrows = "to")
  edges = unique(edges)

  nodes_enriched = unique(node_df[,c(ke_id,"Ke_description",numerical_variables,pval_variable,gene_variable)])
  colnames(nodes_enriched) = c("id","label",numerical_variables,"padj","genes")
  nodes_enriched = unique(nodes_enriched)
  nodes = nodes_enriched

  other_k_idx = which(igraph::V(selegoG)$name %in% nodes_enriched[,"id"] == FALSE) # KE nodes that are not enriched
  if(length(other_k_idx)>0){
    other_nodes = cbind(igraph::V(selegoG)$name[other_k_idx],igraph::V(selegoG)$ke_description[other_k_idx],matrix(NA, nrow = length(other_k_idx),ncol = length(numerical_variables)+2))
    colnames(other_nodes) = c("id","label",numerical_variables,"padj","genes")
    other_nodes = unique(other_nodes)
    nodes = unique(rbind(nodes_enriched, other_nodes))
  }


  for(pi in numerical_variables){
    nodes[[pi]] = as.numeric(as.vector(nodes[[pi]]))
  }

  nodes$padj = as.numeric(as.vector(nodes$padj))

  if(convert_to_gene_symbols){

    # convert node_genes to symbols
    for(i in 1:nrow(nodes)){
      if(length(nodes$genes[i])>0){
        ggi = unlist(strsplit(x = nodes$genes[i], split = ";"))
        ggi_c = convert_genes_to_symbol(genes = ggi,
                                        gene_id_type,
                                        organism,
                                        genes_human=genes_human,
                                        genes_mouse=genes_mouse,
                                        genes_rat=genes_rat)
        nodes$genes[i] = paste(ggi_c,collapse = ";")
      }
    }
  }



  return(list("nodes" = nodes,"edges" = edges))

}

