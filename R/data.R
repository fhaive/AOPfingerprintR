### HUMANS

#' @title genes_human
#' @description Mapping between ENSEMBL, SYMBOL and ENTREZ human gene IDs
#' @format A dataframe with columns ensembl_gene_id, hgnc_symbol, entrezgene_id and description. From human = useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
"genes_human"

#' @title AOP_HUMAN_ENSEMBL_GENES
#' @description Mapping between AOP ids and human ensembl genes
#' @format A list with AOP ids used as indices. Each position contains a vector of human ensembl gene ids
#' @source \url{https://doi.org/10.1038/s41597-023-02321-w}
"human_ens_aop"

#' @title KE_HUMAN_ENSEMBL_GENES
#' @description Mapping between KE ids and human ensembl genes.
#' @format A list with KE ids used as indices. Each position contains a vector of human ensembl gene ids
#' @source \url{https://doi.org/10.1038/s41597-023-02321-w}
"human_ens_ke"

#' @title CLUSTERED_KE_HUMAN_ENSEMBL_GENES
#' @description Mapping between KE ids and human ensembl genes. KEs sharing the same set of genes are clustered in a combined id as KE_id1;KE_id2
#' @format A list with clustered KE ids used as indices. Each position contains a vector of human ensembl gene ids
#' @source \url{https://doi.org/10.1038/s41597-023-02321-w}
"human_ens_clusters"

#' @title AOP_HUMAN_GENE_SYMBOLS
#' @description Mapping between AOP ids and human genes symbols
#' @format A list with AOP ids used as indices. Each position contains a vector of human gene symbol ids
#' @source \url{https://doi.org/10.1038/s41597-023-02321-w}
"human_symbol_aop"

#' @title KE_HUMAN_GENE_SYMBOLS
#' @description Mapping between KE ids and human gene symbols.
#' @format A list with KE ids used as indices. Each position contains a vector of human gene symbol ids
#' @source \url{https://doi.org/10.1038/s41597-023-02321-w}
"human_symbol_ke"

#' @title CLUSTERED_KE_HUMAN_GENE_SYMBOLS
#' @description Mapping between KE ids and human gene symbols. KEs sharing the same set of genes are clustered in a combined id as KE_id1;KE_id2
#' @format A list with clustered KE ids used as indices. Each position contains a vector of human gene symbol ids
#' @source \url{https://doi.org/10.1038/s41597-023-02321-w}
"human_symbol_clusters"

#' @title AOP_HUMAN_ENTREZ_GENES
#' @description Mapping between AOP ids and human genes symbols
#' @format A list with AOP ids used as indices. Each position contains a vector of human gene symbol ids
#' @source \url{https://doi.org/10.1038/s41597-023-02321-w}
"human_entrez_aop"

#' @title KE_HUMAN_ENTREZ_GENES
#' @description Mapping between KE ids and human gene symbols.
#' @format A list with KE ids used as indices. Each position contains a vector of human gene symbol ids
#' @source \url{https://doi.org/10.1038/s41597-023-02321-w}
"human_entrez_ke"

#' @title CLUSTERED_KE_HUMAN_ENTREZ_GENES
#' @description Mapping between KE ids and human gene symbols. KEs sharing the same set of genes are clustered in a combined id as KE_id1;KE_id2
#' @format A list with clustered KE ids used as indices. Each position contains a vector of human gene symbol ids
#' @source \url{https://doi.org/10.1038/s41597-023-02321-w}
"human_entrez_clusters"

#### MOUSE

#' @title genes_mouse
#' @description Mapping between ENSEMBL, SYMBOL and ENTREZ mouse gene IDs
#' @format A dataframe with columns ensembl_gene_id, hgnc_symbol, entrezgene_id and description. From mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
"genes_mouse"

#' @title AOP_MOUSE_ENSEMBL_GENES
#' @description Mapping between AOP ids and mouse ensembl genes
#' @format A list with AOP ids used as indices. Each position contains a vector of mouse ensembl gene ids
#' @source \url{https://doi.org/10.1038/s41597-023-02321-w}
"mouse_ens_aop"

#' @title KE_MOUSE_ENSEMBL_GENES
#' @description Mapping between KE ids and mouse ensembl genes.
#' @format A list with KE ids used as indices. Each position contains a vector of mouse ensembl gene ids
#' @source \url{https://doi.org/10.1038/s41597-023-02321-w}
"mouse_ens_ke"

#' @title CLUSTERED_KE_MOUSE_ENSEMBL_GENES
#' @description Mapping between KE ids and mouse ensembl genes. KEs sharing the same set of genes are clustered in a combined id as KE_id1;KE_id2
#' @format A list with clustered KE ids used as indices. Each position contains a vector of mouse ensembl gene ids
#' @source \url{https://doi.org/10.1038/s41597-023-02321-w}
"mouse_ens_clusters"

#' @title AOP_MOUSE_GENE_SYMBOLS
#' @description Mapping between AOP ids and mouse genes symbols
#' @format A list with AOP ids used as indices. Each position contains a vector of mouse gene symbol ids
#' @source \url{https://doi.org/10.1038/s41597-023-02321-w}
"mouse_symbol_aop"

#' @title KE_MOUSE_GENE_SYMBOLS
#' @description Mapping between KE ids and mouse gene symbols.
#' @format A list with KE ids used as indices. Each position contains a vector of mouse gene symbol ids
#' @source \url{https://doi.org/10.1038/s41597-023-02321-w}
"mouse_symbol_ke"

#' @title CLUSTERED_KE_MOUSE_GENE_SYMBOLS
#' @description Mapping between KE ids and mouse gene symbols. KEs sharing the same set of genes are clustered in a combined id as KE_id1;KE_id2
#' @format A list with clustered KE ids used as indices. Each position contains a vector of mouse gene symbol ids
#' @source \url{https://doi.org/10.1038/s41597-023-02321-w}
"mouse_symbol_clusters"

#' @title AOP_MOUSE_ENTREZ_GENES
#' @description Mapping between AOP ids and mouse genes symbols
#' @format A list with AOP ids used as indices. Each position contains a vector of mouse gene symbol ids
#' @source \url{https://doi.org/10.1038/s41597-023-02321-w}
"mouse_entrez_aop"

#' @title KE_MOUSE_ENTREZ_GENES
#' @description Mapping between KE ids and mouse gene symbols.
#' @format A list with KE ids used as indices. Each position contains a vector of mouse gene symbol ids
#' @source \url{https://doi.org/10.1038/s41597-023-02321-w}
"mouse_entrez_ke"

#' @title CLUSTERED_KE_MOUSE_ENTREZ_GENES
#' @description Mapping between KE ids and mouse gene symbols. KEs sharing the same set of genes are clustered in a combined id as KE_id1;KE_id2
#' @format A list with clustered KE ids used as indices. Each position contains a vector of mouse gene symbol ids
#' @source \url{https://doi.org/10.1038/s41597-023-02321-w}
"mouse_entrez_clusters"

### RAT

#' @title genes_rat
#' @description Mapping between ENSEMBL, SYMBOL and ENTREZ rat gene IDs
#' @format A dataframe with columns ensembl_gene_id, hgnc_symbol, entrezgene_id and description.From  rat = useMart("ensembl", dataset = "rnorvegicus_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
"genes_rat"

#' @title AOP_RAT_ENSEMBL_GENES
#' @description Mapping between AOP ids and rat ensembl genes
#' @format A list with AOP ids used as indices. Each position contains a vector of rat ensembl gene ids
#' @source \url{https://doi.org/10.1038/s41597-023-02321-w}
"rat_ens_aop"

#' @title KE_RAT_ENSEMBL_GENES
#' @description Mapping between KE ids and rat ensembl genes.
#' @format A list with KE ids used as indices. Each position contains a vector of rat ensembl gene ids
#' @source \url{https://doi.org/10.1038/s41597-023-02321-w}
"rat_ens_ke"

#' @title CLUSTERED_KE_RAT_ENSEMBL_GENES
#' @description Mapping between KE ids and rat ensembl genes. KEs sharing the same set of genes are clustered in a combined id as KE_id1;KE_id2
#' @format A list with clustered KE ids used as indices. Each position contains a vector of rat ensembl gene ids
#' @source \url{https://doi.org/10.1038/s41597-023-02321-w}
"rat_ens_clusters"

#' @title AOP_RAT_GENE_SYMBOLS
#' @description Mapping between AOP ids and rat genes symbols
#' @format A list with AOP ids used as indices. Each position contains a vector of rat gene symbol ids
#' @source \url{https://doi.org/10.1038/s41597-023-02321-w}
"rat_symbol_aop"

#' @title KE_RAT_GENE_SYMBOLS
#' @description Mapping between KE ids and rat gene symbols.
#' @format A list with KE ids used as indices. Each position contains a vector of rat gene symbol ids
#' @source \url{https://doi.org/10.1038/s41597-023-02321-w}
"rat_symbol_ke"

#' @title CLUSTERED_KE_RAT_GENE_SYMBOLS
#' @description Mapping between KE ids and rat gene symbols. KEs sharing the same set of genes are clustered in a combined id as KE_id1;KE_id2
#' @format A list with clustered KE ids used as indices. Each position contains a vector of rat gene symbol ids
#' @source \url{https://doi.org/10.1038/s41597-023-02321-w}
"rat_symbol_clusters"

#' @title AOP_RAT_ENTREZ_GENES
#' @description Mapping between AOP ids and rat genes symbols
#' @format A list with AOP ids used as indices. Each position contains a vector of rat gene symbol ids
#' @source \url{https://doi.org/10.1038/s41597-023-02321-w}
"rat_entrez_aop"

#' @title KE_RAT_ENTREZ_GENES
#' @description Mapping between KE ids and rat gene symbols.
#' @format A list with KE ids used as indices. Each position contains a vector of rat gene symbol ids
#' @source \url{https://doi.org/10.1038/s41597-023-02321-w}
"rat_entrez_ke"

#' @title CLUSTERED_KE_RAT_ENTREZ_GENES
#' @description Mapping between KE ids and rat gene symbols. KEs sharing the same set of genes are clustered in a combined id as KE_id1;KE_id2
#' @format A list with clustered KE ids used as indices. Each position contains a vector of rat gene symbol ids
#' @source \url{https://doi.org/10.1038/s41597-023-02321-w}
"rat_entrez_clusters"


### AOPs
#' @title AOP_NAMES
#' @description AOP_IDs_NAMES_ASSOCIATION
#' @format A data frame with 2 variables:
#' \describe{
#'   \item{\code{a.AOP_ID}}{AOP ids}
#'   \item{\code{a.name}}{AOP names}
#'}
#' @source \url{https://doi.org/10.1038/s41597-023-02321-w}
"aop_names"

#' @title BIOLOGICAL_SYSTEM
#' @description BIOLOGICAL_SYSTEM_ANNOTATION
#' @format A data frame with 11 variables:
#' \describe{
#'   \item{\code{a.AOP_ID}}{AOP ids}
#'   \item{\code{a.name}}{AOP names}
#'   \item{\code{key_event_name}}{Name of the key event in the biological system}
#'   \item{\code{ke}}{Key event identifier}
#'   \item{\code{level}}{Biological level of organization}
#'   \item{\code{system}}{Primary biological system involved}
#'   \item{\code{organ_tissue}}{Primary organ or tissue involved}
#'   \item{\code{cell}}{Primary cell type involved}
#'   \item{\code{cell_component}}{Primary cell component involved}
#'   \item{\code{secondary_system}}{Secondary biological system involved, if any}
#'   \item{\code{secondary_organ_tissue}}{Secondary organ or tissue involved, if any}
#'   \item{\code{secondary_cell}}{Secondary cell type involved, if any}
#'   \item{\code{secondary_cell_component}}{Secondary cell component involved, if any}
#'}
#' @source \url{https://doi.org/10.1038/s41597-023-02321-w}
"Biological_system_annotations"

### AOPs SSBD ANNOTATION
#' @title AOP_SSbD_ANNOTATION
#' @description AOP_SSbD_IMPACT_ASSOCIATION
#' @format A data frame with 6 variables:
#' \describe{
#'   \item{\code{AOP}}{AOP ids}
#'   \item{\code{SSbD_category}}{SSbD impact category}
#'   \item{\code{Endpoint}}{Endpoint}
#'   \item{\code{Organ}}{Organ}
#'   \item{\code{AOP_name}}{AOP_name}
#'   \item{\code{Notes}}{Notes}
#'}
#' @source \url{https://doi.org/10.1038/s41597-023-02321-w}
"Annotate_AOPs"

#' @title KE_NAMES
#' @description KE_IDs_NAMES_ASSOCIATION
#' @format A data frame with 2 variables:
#' \describe{
#'   \item{\code{Ke}}{Ke ids}
#'   \item{\code{Ke_description}}{KE names}
#'}
#' @source \url{https://doi.org/10.1038/s41597-023-02321-w}
"ke_names"

#' @title CLUSTERED_KE
#' @description KE-cluster associations
#' @format A data frame with 3 variables:
#' \describe{
#'   \item{\code{cluster}}{KE cluster ids}
#'   \item{\code{KEs}}{KE ids}
#'   \item{\code{Description}}{KE names}
#'}
#' @source \url{https://doi.org/10.1038/s41597-023-02321-w}
"clust_anno"

#' @title AOP_KE_ASSOCIATION
#' @description AOP_KE_MAPPING
#' @format A data frame with 5 variables:
#' \describe{
#'   \item{\code{Aop}}{AOP ids}
#'   \item{\code{Ke}}{KE ids}
#'   \item{\code{Ke_type}}{KE type description, possible values are MolecularInitiatingEvent, AdverseOutcome, KeyEvent}
#'   \item{\code{Ke_description}}{KE name}
#'   \item{\code{a.name}}{AOP name}
#'}
#' @source \url{https://doi.org/10.1038/s41597-023-02321-w}
"aop_ke_table_hure"

#' @title AOP_KE_GENE_ASSOCIATION_DF
#' @description Dataframe of association between AOPs, KEs and ensembl genes
#' @format A data frame with 6 variables:
#' \describe{
#'   \item{\code{Aop}}{AOP ids}
#'   \item{\code{Aop_KE}}{KE ids belonging to Aop id}
#'   \item{\code{KE}}{KE id}
#'   \item{\code{Annotation}}{GO id mapped to KE}
#'   \item{\code{Level}}{Level}
#'   \item{\code{Gene}}{Gene id}
#'}
#' @source \url{https://doi.org/10.1038/s41597-023-02321-w}
"aops_new"

#' @title KE_KE_NETWORK
#' @description KE-KE network.
#' @format It is an igraph object, whose nodes are the clustered KEs, and edges represent their connections. Nodes have the attribute 'ke_description' that includes a description of the KEs, while edges have the attribute 'type.r' that specifies the type of connections.
#' @source \url{https://doi.org/10.1038/s41597-023-02321-w}
"KEKE_net"
