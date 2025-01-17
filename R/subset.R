#' Subset Herbs and Retrieve Associated Molecules and Targets
#'
#' This function retrieves herb, molecule, and target information from the \code{Yulab_data} dataset based on specified herb names and their corresponding name types (Chinese, Pinyin, or English).
#'
#' @param herb A character vector containing the names of herbs to be queried.
#' @param type A string indicating the type of herb names provided. Must be one of \code{"Herb_cn_name"}, \code{"Herb_pinyin_name"}, or \code{"Herb_en_name"}.
#'
#' @return A \code{data.frame} with three columns: \code{herb}, \code{molecule}, and \code{target}, containing the corresponding information for the specified herbs. Rows with \code{NA} values are excluded.
#'
#' @import dplyr
#' @import tidyr
#' @examples
#' \dontrun{
#' # Example usage with Chinese herb names
#' herbs <- c("灵芝")
#' lz <- subset_herb(herb = herbs, type = "Herb_cn_name")
#' print(lz)
#'
#' # Example usage with pinyin herb names
#' herbs <- c("Ginseng", "Licorice")
#' comp <- subset_herb(herb = herbs, type = "Herb_pinyin_name")
#' print(comp)
#' }
#'
#' @export
#' 
subset_herb <- function(herb, type){
  
  # Load necessary packages
  library(dplyr)
  library(tidyr)
  load("Yulab_data.rda")
  
  # Validate the 'type' parameter
  type <- match.arg(type, c("Herb_cn_name", "Herb_pinyin_name", "Herb_en_name"))
  
  # Check existence of herbs based on the specified type
  if (type == "Herb_cn_name"){
    if (all(herb %in% unique(Yulab_data$Herb_cn_name)) == FALSE){
      herb_not_exist <- setdiff(herb, unique(Yulab_data$Herb_cn_name))
      print(paste0(herb_not_exist, " doesn't/don't exist in our dataset."))
      herb <- herb[-match(herb_not_exist, herb)]
    }
    
    result <- Yulab_data %>%
      dplyr::filter(Herb_cn_name %in% herb) %>%
      dplyr::select(c("Herb_pinyin_name", "molecule", "target")) %>%
      dplyr::distinct(.keep_all = TRUE) %>%
      tidyr::drop_na()
    colnames(result) <- c("herb", "molecule", "target")
    return(result)
  }
  
  else if (type == "Herb_pinyin_name"){
    if (all(herb %in% unique(Yulab_data$Herb_pinyin_name)) == FALSE){
      herb_not_exist <- setdiff(herb, unique(Yulab_data$Herb_pinyin_name))
      print(paste0(herb_not_exist, " doesn't/don't exist in our dataset."))
      herb <- herb[-match(herb_not_exist, herb)]
    }
    
    result <- Yulab_data %>%
      dplyr::filter(Herb_pinyin_name %in% herb) %>%
      dplyr::select(c("Herb_pinyin_name", "molecule", "target")) %>%
      dplyr::distinct(.keep_all = TRUE) %>%
      tidyr::drop_na()
    colnames(result) <- c("herb", "molecule", "target")
    return(result)
  }
  
  else if (type == "Herb_en_name"){
    if (all(herb %in% unique(Yulab_data$Herb_en_name)) == FALSE){
      herb_not_exist <- setdiff(herb, unique(Yulab_data$Herb_en_name))
      print(paste0(herb_not_exist, " doesn't/don't exist in our dataset."))
      herb <- herb[-match(herb_not_exist, herb)]
    }
    
    result <- Yulab_data %>%
      dplyr::filter(Herb_en_name %in% herb) %>%
      dplyr::select("Herb_pinyin_name", "molecule", "target") %>%
      dplyr::distinct(.keep_all = TRUE) %>%
      tidyr::drop_na()
    colnames(result) <- c("herb", "molecule", "target")
    return(result)
  }
  
  else{
    stop("Parameter 'type' should be one of 'Herb_cn_name', 'Herb_pinyin_name', or 'Herb_en_name'.")
  }
}


#' Retrieve Herbs, Molecules, and Targets Based on Gene List
#'
#' This function retrieves herb, molecule, and target information from the \code{Yulab_data} dataset based on a provided list of target genes.
#'
#' @param gene_list A character vector containing gene symbols which can be determined freely by users.
#'
#' @return A \code{data.frame} with three columns: \code{herb}, \code{molecule}, and \code{target}, containing information related to the specified genes. Rows with \code{NA} values are excluded.
#'
#' @import dplyr
#' @import tidyr
#' @examples
#' \dontrun{
#' # Example usage with a list of gene symbols
#' genes <- c("TP53", "EGFR", "BRCA1")
#' herbs_targets <- target_search(genes)
#' print(herbs_targets)
#' }
#'
#' @export
#' 
target_search <- function(gene_list){
  
  # Load necessary packages
  library(dplyr)
  library(tidydr)
  load("Yulab_data.rda")
  
  # Check existence of genes in the dataset
  if (!all(gene_list %in% unique(Yulab_data$target))){
    gene_diff <- setdiff(gene_list, unique(Yulab_data$target))
    print(paste0(gene_diff, " doesn't/don't exist in the datasets."))
    gene_list <- gene_list[-match(gene_diff, gene_list)]
  }
  
  # Retrieve relevant data
  herbs_data <- Yulab_data %>%
    dplyr::filter(target %in% gene_list) %>%
    dplyr::select(c("Herb_pinyin_name", "molecule", "target")) %>%
    dplyr::distinct(.keep_all = TRUE) %>%
    tidyr::drop_na()
  
  colnames(herbs_data) <- c("herb", "molecule", "target")
  return(herbs_data)
}
