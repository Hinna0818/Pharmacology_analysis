#' Perform KEGG Pathway Enrichment Analysis on Target Genes
#'
#' This function performs KEGG pathway enrichment analysis using the \code{clusterProfiler} package on a set of target genes. It takes a data frame containing target genes and optionally a specific gene list to analyze.
#'
#' @param x A \code{data.frame} containing at least a \code{target} column with gene symbols.
#' @param gene_list A character vector of gene symbols to be used for enrichment analysis. Defaults to all targets in \code{x}.
#' @param organism A string specifying the organism code for KEGG analysis. Default is \code{"hsa"} for Homo sapiens.
#' @param OrgDb A string specifying the OrgDb object for gene ID conversion. Default is \code{"org.Hs.eg.db"}.
#' @param pvalueCutoff A numeric value indicating the p-value cutoff for significance. Default is \code{0.05}.
#' @param pAdjustMethod A string specifying the multiple testing correction method. Default is \code{"BH"} (Benjamini-Hochberg).
#' @param qvalueCutoff A numeric value indicating the q-value cutoff for significance. Default is \code{0.2}.
#'
#' @return An \code{enrichResult} object containing the KEGG enrichment analysis results.
#'
#' @import clusterProfiler
#' @import dplyr
#' @import OrgDb
#'
#' @examples
#' \dontrun{
#' # using data from "fengmi" as an example
#'fm <- subset_herb("fengmi", type = "Herb_pinyin_name")
#'random <- sample(1:nrow(fm), 100, replace = FALSE)
#'sub_data <- fm[random, ]
#'gene_list <- sub_data$target
#'fm_enrich <- enrich_target_KEGG(x = fm, gene_list = gene_list)
#'clusterProfiler::dotplot(fm_enrich, title = "KEGG")
#'clusterProfiler::cnetplot(fm_enrich)
#' }
#'
#' @export
#' 
enrich_target_KEGG <- function(
    x, ## A data.frame containing target gene information
    gene_list = x$target, ## A gene set for enrichment analysis; defaults to all targets in x
    organism = "hsa",
    OrgDb = "org.Hs.eg.db",
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    qvalueCutoff = 0.2
){
  
  # Load necessary package
  library(clusterProfiler)
  
  # Validate input data
  if (!"target" %in% colnames(x)){
    stop("Input data frame 'x' must contain a 'target' column with gene symbols.")
  }
  
  # Intersect gene lists
  gene_input <- intersect(x$target, gene_list)
  
  if(length(gene_input) == 0){
    stop("No overlapping genes between 'x$target' and 'gene_list'.")
  }
  
  # Convert gene symbols to ENTREZID
  eg <- bitr(unique(gene_input), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = OrgDb)
  
  # Perform KEGG enrichment analysis
  enrich_result <- enrichKEGG(
    gene = eg$ENTREZID,
    organism = organism,
    pvalueCutoff = pvalueCutoff,
    pAdjustMethod = pAdjustMethod,
    qvalueCutoff = qvalueCutoff
  )
  
  # Convert ENTREZID back to readable gene symbols
  enrich_result <- setReadable(enrich_result, OrgDb, keyType = "ENTREZID")
  
  return(enrich_result)
  
}


#' Perform Gene Ontology (GO) Enrichment Analysis on Target Genes
#'
#' This function conducts Gene Ontology (GO) enrichment analysis using the \code{clusterProfiler} package on a set of target genes. It allows specification of the GO ontology category (BP, MF, or CC) and other enrichment parameters.
#'
#' @param x A \code{data.frame} containing at least a \code{target} column with gene symbols.
#' @param gene_list A character vector of gene symbols to be used for enrichment analysis. Defaults to all targets in \code{x}.
#' @param ont A string specifying the GO ontology category. Must be one of \code{"BP"} (Biological Process), \code{"MF"} (Molecular Function), or \code{"CC"} (Cellular Component). Default is \code{"BP"}.
#' @param OrgDb A string specifying the OrgDb object for gene ID conversion. Default is \code{"org.Hs.eg.db"}.
#' @param pvalueCutoff A numeric value indicating the p-value cutoff for significance. Default is \code{0.05}.
#' @param pAdjustMethod A string specifying the multiple testing correction method. Default is \code{"BH"} (Benjamini-Hochberg).
#' @param qvalueCutoff A numeric value indicating the q-value cutoff for significance. Default is \code{0.2}.
#' @param readable A logical indicating whether to convert gene IDs to gene symbols in the output. Default is \code{TRUE}.
#'
#' @return An \code{enrichResult} object containing the GO enrichment analysis results.
#'
#' @import clusterProfiler
#' @import dplyr
#' @import OrgDb
#'
#' @examples
#' \dontrun{
#' # using data from "chantui" for an example
#'ct <- subset_herb(herb = "蝉蜕", type = "Herb_cn_name")
#'ct_enrich <- enrich_target_GO(x = ct, gene_list = ct$target) ## using all targets in chantui as default
#'clusterProfiler::dotplot(ct_enrich,title = "BP")
#' }
#'
#' @export
enrich_target_GO <- function(
    x, ## A data.frame containing target gene information
    gene_list = x$target, ## A gene set for enrichment analysis; defaults to all targets in x
    ont = "BP",
    OrgDb = "org.Hs.eg.db",
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    qvalueCutoff = 0.2,
    readable = TRUE
){
  
  # Load necessary package
  library(clusterProfiler)
  
  # Validate input data
  if (!"target" %in% colnames(x)){
    stop("Input data frame 'x' must contain a 'target' column with gene symbols.")
  }
  
  # Validate 'ont' parameter
  ont <- match.arg(ont, c("BP", "MF", "CC"))
  
  # Intersect gene lists
  gene_input <- intersect(x$target, gene_list)
  
  if(length(gene_input) == 0){
    stop("No overlapping genes between 'x$target' and 'gene_list'.")
  }
  
  # Convert gene symbols to ENTREZID
  eg <- bitr(unique(gene_input), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = OrgDb)
  
  # Perform GO enrichment analysis
  enrich_result <- enrichGO(
    gene = eg$ENTREZID,
    OrgDb = OrgDb,
    ont = ont,
    pvalueCutoff = pvalueCutoff,
    pAdjustMethod = pAdjustMethod,
    qvalueCutoff = qvalueCutoff
  )
  
  # Optionally convert ENTREZID back to readable gene symbols
  if(readable){
    enrich_result <- setReadable(enrich_result, OrgDb, keyType = "ENTREZID")
  }
  
  return(enrich_result)
}


#' Perform Disease Ontology (DO) Enrichment Analysis on Target Genes
#'
#' This function conducts Disease Ontology (DO) enrichment analysis using the \code{DOSE} package on a set of target genes. It takes a data frame containing target genes and optionally a specific gene list to analyze.
#'
#' @param x A \code{data.frame} containing at least a \code{target} column with gene symbols.
#' @param gene_list A character vector of gene symbols to be used for enrichment analysis. Defaults to all targets in \code{x}.
#' @param ont A string specifying the DO ontology category. Must be one of \code{"HDO"} (Human Disease Ontology), \code{"HPO"} (Human Phenotype Ontology), or \code{"MPO"} (Mouse Phenotype Ontology). Default is \code{"HDO"}.
#' @param OrgDb A string specifying the OrgDb object for gene ID conversion. Default is \code{"org.Hs.eg.db"}.
#' @param pvalueCutoff A numeric value indicating the p-value cutoff for significance. Default is \code{0.05}.
#' @param pAdjustMethod A string specifying the multiple testing correction method. Default is \code{"BH"} (Benjamini-Hochberg).
#' @param qvalueCutoff A numeric value indicating the q-value cutoff for significance. Default is \code{0.2}.
#' @param readable A logical indicating whether to convert gene IDs to gene symbols in the output. Default is \code{TRUE}.
#'
#' @return An \code{enrichResult} object containing the Disease Ontology enrichment analysis results.
#'
#' @import DOSE
#' @import clusterProfiler
#' @import dplyr
#' @import org.Hs.eg.db
#' @import enrichplot
#'
#' @examples
#' \dontrun{
#' # Example usage with a data frame containing target genes
#' df <- data.frame(target = c("TP53", "EGFR", "BRCA1", "MYC", "CDK2", "MDR1", "STAT3", "AKT1", "MTOR", "BCL2"))
#' do_enrich_results <- enrich_target_DO(df)
#' print(do_enrich_results)
#'
#' # Visualize the enrichment results
#' barplot(do_enrich_results, showCategory=10, title = "Disease Ontology Enrichment")
#' dotplot(do_enrich_results, showCategory=10, title = "Disease Ontology Enrichment")
#' cnetplot(do_enrich_results, showCategory=5)
#' emapplot(do_enrich_results, showCategory=5)
#' }
#'
#' @export
enrich_target_DO <- function(
    x, ## A data.frame containing target gene information
    gene_list = x$target, ## A gene set for enrichment analysis; defaults to all targets in x
    ont = "HDO", ## using HDO(Human Disease Ontology) as default
    OrgDb = "org.Hs.eg.db",
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    qvalueCutoff = 0.2,
    readable = TRUE
){
  
  # Load necessary packages
  library(DOSE)
  library(clusterProfiler)
  
  # Validate input data
  if (!"target" %in% colnames(x)){
    stop("Input data frame 'x' must contain a 'target' column with gene symbols.")
  }
  
  ont <- match.arg(ont, c("HDO", "HPO", "MPO"))
  
  # Intersect gene lists
  gene_input <- intersect(x$target, gene_list)
  
  if(length(gene_input) == 0){
    stop("No overlapping genes between 'x$target' and 'gene_list'.")
  }
  
  # Convert gene symbols to ENTREZID
  eg <- bitr(unique(gene_input), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = OrgDb)
  
  # Perform Disease Ontology enrichment analysis
  do_enrich_result <- enrichDO(
    gene = eg$ENTREZID,
    ont = "HDO",                
    pvalueCutoff = pvalueCutoff,
    pAdjustMethod = pAdjustMethod,
    minGSSize = 10,
    maxGSSize = 500,
    readable = readable        # Convert ENTREZID to gene symbols if TRUE
  )
  
  return(do_enrich_result)
  
}
