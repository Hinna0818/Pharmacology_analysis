# Pharmacology Analysis 🔬🌿

The `pharmacology_analysis` repository provides an integrated pipeline for Traditional Chinese Medicine (TCM) network pharmacology, including compound-target data mining, functional enrichment analysis, and network visualization.

## 🌟 Features

- 🔍 **TCM Compound–Target Retrieval**  
  Curate relationships between herbs, their bioactive compounds, and predicted/known targets.

- 📊 **Target Enrichment Analysis**  
  Perform GO and KEGG pathway enrichment using `clusterProfiler`.

- 🌐 **Network Visualization**  
  Visualize Herb–Compound–Target networks using `igraph` and `ggraph` in customizable concentric layouts.

## 📦 Dependencies

```r
# Suggested R version: ≥ 4.2
install.packages(c("tidyverse", "ggraph", "igraph", "clusterProfiler", 
                   "org.Hs.eg.db", "RColorBrewer"))
```

## 🧪 Applications
This framework is useful for:
- Exploring pharmacological mechanisms of TCM
- Identifying hub genes or druggable targets
- Integrating multi-omics and network pharmacology
