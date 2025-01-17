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
#' @param minGSSize A numeric value indicating the minimum size of gene sets in enrichment analysis. Default is \code{10}.
#' @param maxGSSize A numeric value indicating the maximum size of gene sets in enrichment analysis. Default is \code{500}.
#'
#' @return An \code{enrichResult} object containing the KEGG enrichment analysis results.
#'
#' @import clusterProfiler
#' @import dplyr
#' @import org.Hs.eg.db
#'
#' @examples
#' \dontrun{
#' # using data from "fengmi" as an example
#'fm <- subset_herb(herb = "fengmi", type = "Herb_pinyin_name")
#'fm_enrich <- enrich_target_KEGG(x = fm, gene_list = fm$target)
#'
## dotplot
#'clusterProfiler::dotplot(fm_enrich, title = "KEGG Pathway Enrichment")
#'
## barplot
#'barplot(fm_enrich, showCategory=10, title = "KEGG Pathway Enrichment") + theme_minimal()
#'
## cnetplot
#'clusterProfiler::cnetplot(fm_enrich, showCategory=3, foldChange=NULL, circular=FALSE, colorEdge=TRUE) + theme_minimal()
#' }
#'
#' @export
#' 
enrich_target_KEGG <- function(
    x, 
    gene_list = x$target, 
    organism = "hsa",
    OrgDb = org.Hs.eg.db,
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    qvalueCutoff = 0.2,
    minGSSize = 10,
    maxGSSize = 500
){
  
  library(clusterProfiler)
  library(org.Hs.eg.db)
  
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
  
  if(nrow(eg) == 0){
    stop("No gene symbols could be converted to ENTREZID.")
  }
  
  # Perform KEGG enrichment analysis
  enrich_result <- enrichKEGG(
    gene = eg$ENTREZID,
    organism = organism,
    pvalueCutoff = pvalueCutoff,
    pAdjustMethod = pAdjustMethod,
    qvalueCutoff = qvalueCutoff,
    minGSSize = minGSSize,
    maxGSSize = maxGSSize
  )
  
  # Check if enrichment result is not empty
  if(is.null(enrich_result) || nrow(enrich_result@result) == 0){
    warning("No significant KEGG pathways found.")
    return(enrich_result)
  }
  
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
#' @param minGSSize A numeric value indicating the minimum size of gene sets in enrichment analysis. Default is \code{10}.
#' @param maxGSSize A numeric value indicating the maximum size of gene sets in enrichment analysis. Default is \code{500}.
#' @param readable A logical indicating whether to convert gene IDs to gene symbols in the output. Default is \code{TRUE}.
#'
#' @return An \code{enrichResult} object containing the GO enrichment analysis results.
#'
#' @import clusterProfiler
#' @import dplyr
#' @import org.Hs.eg.db
#'
#' @examples
#' \dontrun{
#' using data from "chantui" for an example
#'ct <- subset_herb(herb = "蝉蜕", type = "Herb_cn_name")
#'ct_enrich <- enrich_target_GO(x = ct, gene_list = ct$target, ont = "BP")
#'

## dotplot
#'clusterProfiler::dotplot(ct_enrich, title = "GO Pathway Enrichment")
#'

## barplot
#'barplot(ct_enrich, showCategory=10, title = "GO Pathway Enrichment") + theme_minimal()
#'

## cnetplot
#'clusterProfiler::cnetplot(ct_enrich, showCategory=3, foldChange=NULL, circular=FALSE, colorEdge=TRUE) + theme_minimal()
#' }
#'
#' @export
#' 
enrich_target_GO <- function(
    x, 
    gene_list = x$target, 
    ont = "BP",
    OrgDb = org.Hs.eg.db,
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    qvalueCutoff = 0.2,
    minGSSize = 10,
    maxGSSize = 500,
    readable = TRUE
){
  
  library(clusterProfiler)
  library(org.Hs.eg.db)
  
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
  eg <- clusterProfiler::bitr(unique(gene_input), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = OrgDb)
  
  if(nrow(eg) == 0){
    stop("No gene symbols could be converted to ENTREZID.")
  }
  
  # Perform GO enrichment analysis
  enrich_result <- clusterProfiler::enrichGO(
    gene = eg$ENTREZID,
    OrgDb = OrgDb,
    ont = ont,
    pvalueCutoff = pvalueCutoff,
    pAdjustMethod = pAdjustMethod,
    qvalueCutoff = qvalueCutoff,
    minGSSize = minGSSize,
    maxGSSize = maxGSSize,
    readable = readable
  )
  
  # Check if enrichment result is not empty
  if(is.null(enrich_result) || nrow(enrich_result@result) == 0){
    warning("No significant GO terms found.")
    return(enrich_result)
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
#' @param organism A string specifying the organism type. Default is \code{"hsa"}.
#' @param pvalueCutoff A numeric value indicating the p-value cutoff for significance. Default is \code{0.05}.
#' @param pAdjustMethod A string specifying the multiple testing correction method. Default is \code{"BH"} (Benjamini-Hochberg).
#' @param minGSSize A numeric value indicating the minimum size of gene sets in enrichment analysis. Default is \code{10}.
#' @param maxGSSize A numeric value indicating the maximum size of gene sets in enrichment analysis. Default is \code{500}.
#' @param qvalueCutoff A numeric value indicating the q-value cutoff for significance. Default is \code{0.2}.
#' @param readable A logical indicating whether to convert gene IDs to gene symbols in the output. Default is \code{TRUE}.
#'
#' @return An \code{enrichResult} object containing the Disease Ontology enrichment analysis results.
#'
#' @import DOSE
#' @import clusterProfiler
#' @import dplyr
#' @import org.Hs.eg.db
#'
#' @examples
#' \dontrun{
#' # using "chantui" as example
#'ct_enrich_DO <- enrich_target_DO(x = ct, ont = "HDO", gene_list = ct$target)
#'

## dotplot
#'clusterProfiler::dotplot(ct_enrich_DO, title = "DO Pathway Enrichment")
#'

## barplot
#'barplot(ct_enrich_DO, showCategory=10, title = "DO Pathway Enrichment") + theme_minimal()
#'

## cnetplot
#'clusterProfiler::cnetplot(ct_enrich_DO, showCategory=3, foldChange=NULL, circular=FALSE, colorEdge=TRUE) + theme_minimal()
#' }
#'
#' @export
#' 
enrich_target_DO <- function(
    x, 
    gene_list = x$target, 
    ont = "HDO", ## using HDO(Human Disease Ontology) as default
    OrgDb = org.Hs.eg.db::org.Hs.eg.db,
    organism = "hsa",
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    minGSSize = 10,
    maxGSSize = 500,
    qvalueCutoff = 0.2,
    readable = TRUE
){
  
  library(clusterProfiler)
  library(org.Hs.eg.db)
  
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
  eg <- clusterProfiler::bitr(unique(gene_input), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = OrgDb)
  
  if(nrow(eg) == 0){
    stop("No gene symbols could be converted to ENTREZID.")
  }
  
  # Perform Disease Ontology enrichment analysis
  do_enrich_result <- DOSE::enrichDO(
    gene = eg$ENTREZID,
    ont = ont,
    organism = organism,
    pvalueCutoff = pvalueCutoff,
    pAdjustMethod = pAdjustMethod,
    qvalueCutoff = qvalueCutoff,
    minGSSize = 10,
    maxGSSize = 500,
    readable = readable        
  )
  
  # Check if enrichment result is not empty
  if(is.null(do_enrich_result) || nrow(do_enrich_result@result) == 0){
    warning("No significant Disease Ontology terms found.")
    return(do_enrich_result)
  }
  
  return(do_enrich_result)
  
}
