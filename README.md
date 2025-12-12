# AOPfingerprintR

**AOPfingerprintR** is an R package that provides curated lists of gene identifiers mapped to **Adverse Outcome Pathways (AOPs)** and their **Key Events (KEs)**, as described in [DOI: 10.1038/s41597-023-02321-w](https://doi.org/10.1038/s41597-023-02321-w).

## Features

- **Gene identifiers**: Supports multiple identifier types, including  
  - Ensembl  
  - Gene Symbol  
  - Entrez  

- **Enrichment functions**: Tools to enrich **Key Events (KEs)** and **Adverse Outcome Pathways (AOPs)**.  

- **AOPfingerprint analysis**: Implements the methodology described in  
  [DOI: 10.1002/advs.202203984](https://doi.org/10.1002/advs.202203984).  

- **Visualization**: Functions to plot **KE–KE interaction networks** resulting from KE enrichment analysis.  

---

## References

- Data source: [10.1038/s41597-023-02321-w](https://doi.org/10.1038/s41597-023-02321-w)  
- Methodology: [10.1002/advs.202203984](https://doi.org/10.1002/advs.202203984)

---

## 📖 Manual
The package manual is available here:  
[**AOPfingerprintR_0.0.1.pdf**](https://github.com/fhaive/AOPfingerprintR/blob/main/AOPfingerprintR_0.0.1.pdf)

The manual was generated using:
```r
devtools::build_manual(path = ".")
```

---

## 💻 Installation (from GitHub)

1. **Get "remotes" package**:

   ```r
   install.packages("remotes")
   ```
   
2. **Install the package from GitHub**
   
   ```r
   remotes::install_github(repo="fhaive/AOPfingerprintR")
   ```
   
## 💻 Installation (build local)

1. **Create tar.gz file**:  

   ```r
   devtools::build()
   ```

2. **Install the package**:
   ```r
   install.packages("path_to_the_tar.gz.file", type = "source", repos = NULL)
   ```
   
---

## 👩‍🔬 Authors & Acknowledgments
- **Angela Serra**  
- **Michele Fratello**  

---

## 📌 Project Status
The project is **ongoing**.

---   

## 📊 Sample Data
Example data for bleomycin exposure at multiple doses and time points is available here:  
[**Bleomycin Data**](https://github.com/fhaive/AOPfingerprintR/tree/main/sample_inputs)

## 📝 Example Usage of *AOPfingerprintR*

The following examples demonstrate how to use **AOPfingerprintR** depending on whether your input gene lists include numerical variables (e.g., BMD values) or not.

---

### 1️⃣ Example without numerical variables

In this case, the uploaded gene lists contain only gene identifiers, with no numerical variables.

```r
library(AOPfingerprintR)

# Load background genes
background <- readxl::read_excel("sample_inputs/background_genes.xlsx", col_names = FALSE)
background <- background$...1

# Load gene lists from Excel
file <- "sample_inputs/dose_dependent_genes_bleomycin_no_PODs.xlsx"
sheets <- readxl::excel_sheets(file)
GList <- lapply(sheets, function(X) as.data.frame(readxl::read_excel(file, sheet = X)))
names(GList) <- sheets

# KE enrichment
numerical_properties <- NULL
ke_enrichment_results <- enrich_KEs_AOPs(
  GList = GList,
  list_gene_sets = human_ens_clusters,
  only_significant = TRUE,
  pval_th = 0.05,
  adj.method = "fdr",
  merge_by = "Ke",
  numerical_properties = numerical_properties,
  background = background
)

# AOP enrichment
aop_enrichment_results <- enrich_KEs_AOPs(
  GList = GList,
  list_gene_sets = human_ens_aop,
  only_significant = TRUE,
  pval_th = 0.05,
  adj.method = "fdr",
  merge_by = "Aop",
  numerical_properties = numerical_properties,
  background = background
)

# Build AOP fingerprints
res <- build_aop_for_aop_fingeprints(
  aop_enrichment_results,
  ke_enrichment_results,
  min_aop_length = 5,
  percentage_enriched_ke = 0.33
)

# Bubble plot visualization
library(ggplot2)
p <- render_aop_fingerprint_bubble_plot(
  enrichement_data = res$detailed_results_only_enriched,
  x_axis_var = "Experiment",
  y_axis_var = "AOP_name",
  threshold_proportion = 0.33,
  text_cex = 15
)
p

# KE–KE interaction network visualization
detailed_results <- ke_enrichment_results
nodes_edges <- make_visNetwork(
  detailed_results,
  experiment = "bleomycin_72",
  enlarge_ke_selection = TRUE,
  ke_id = "TermID",
  numerical_variables = numerical_properties,
  pval_variable = "padj",
  gene_variable = "Genes",
  max_path_length = 3,
  n_AOs = 2,
  n_MIEs = 2
)

vn <- plot_visNetwork(
  nodes = nodes_edges$nodes,
  edges = nodes_edges$edges,
  group_by = "ssbd",
  numerical_variables = numerical_properties
)

vn
```


### 2️⃣ Example with numerical variables

In this case, each gene in the input lists has associated BMD, BMDL, and BMDU values.

```r
library(AOPfingerprintR)

# Load background genes
background <- readxl::read_excel("sample_inputs/background_genes.xlsx", col_names = FALSE)
background <- background$...1

# Load gene lists from Excel
file <- "sample_inputs/dose_dependent_genes_bleomycin.xlsx"
sheets <- readxl::excel_sheets(file)
GList <- lapply(sheets, function(X) as.data.frame(readxl::read_excel(file, sheet = X)))
names(GList) <- sheets

# KE enrichment with numerical properties
ke_enrichment_results <- enrich_KEs_AOPs(
  GList = GList,
  list_gene_sets = human_ens_clusters,
  only_significant = TRUE,
  pval_th = 0.05,
  adj.method = "fdr",
  merge_by = "Ke",
  numerical_properties = c("BMD","BMDL","BMDU"),
  background = background
)

# AOP enrichment with numerical properties
aop_enrichment_results <- enrich_KEs_AOPs(
  GList = GList,
  list_gene_sets = human_ens_aop,
  only_significant = TRUE,
  pval_th = 0.05,
  adj.method = "fdr",
  merge_by = "Aop",
  numerical_properties = c("BMD","BMDL","BMDU"),
  background = background
)

# Build AOP fingerprints
res <- build_aop_for_aop_fingeprints(
  aop_enrichment_results,
  ke_enrichment_results,
  min_aop_length = 5,
  percentage_enriched_ke = 0.33
)

# Bubble plot visualization
library(ggplot2)
p <- render_aop_fingerprint_bubble_plot(
  enrichement_data = res$detailed_results_only_enriched,
  x_axis_var = "Experiment",
  y_axis_var = "AOP_name",
  filter_column = NULL,
  filter_by = "STO_tox_Lung",
  threshold_proportion = 0.33,
  text_cex = 15
)
p

# KE–KE interaction network visualization
detailed_results <- ke_enrichment_results
nodes_edges <- make_visNetwork(
  detailed_results,
  experiment = "bleomycin_72",
  enlarge_ke_selection = TRUE,
  ke_id = "TermID",
  numerical_variables = c("BMDL","BMD","BMDU"),
  pval_variable = "padj",
  gene_variable = "Genes",
  max_path_length = 3,
  n_AOs = 2,
  n_MIEs = 2
)

vn <- plot_visNetwork(
  nodes = nodes_edges$nodes,
  edges = nodes_edges$edges,
  group_by = "ssbd",
  numerical_variables = c("BMDL","BMD","BMDU")
)
vn

```

