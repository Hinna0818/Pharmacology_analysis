## 根据特定中药名查找对应单体和靶点信息函数
subset_herb <- function(herb, type){
  
  # 加载数据
  library(dplyr)
  load("Yulab_data.rda")
  
  # 判断type的值，根据拼音还是中文还是英文进行查询
  type <- match.arg(type, c("Herb_cn_name", "Herb_pinyin_name", "Herb_en_name"))
  
  # 如果某个中药在数据库中不存在，那就只返回存在的查询结果
  if (type == "Herb_cn_name"){
    if (all(herb %in% (unique(Yulab_data$Herb_cn_name))) == FALSE){
      herb_not_exist <- setdiff(herb, unique(Yulab_data$Herb_cn_name))
      print(paste0(herb_not_exist, "doesn't/don't exist in our Yulab_datas."))
      herb <- herb[-match(herb_not_exist, herb)]
    }
    
    result <- Yulab_data[which((Yulab_data$Herb_cn_name) %in% herb), ] %>%
      dplyr::select(c(2, 4, 5)) %>%
      dplyr::distinct(.keep_all = TRUE) %>%
      tidyr::drop_na()
    colnames(result) <- c("herb", "molecule", "target")
    return(result)
  }
  
  else if (type == "Herb_pinyin_name"){
    if (all(herb %in% (unique(Yulab_data$Herb_pinyin_name))) == FALSE){
      herb_not_exist <- setdiff(herb, unique(Yulab_data$Herb_pinyin_name))
      print(paste0(herb_not_exist, "doesn't/don't exist in our Yulab_datas."))
      herb <- herb[-match(herb_not_exist, herb)]
    }
    
    result <- Yulab_data[which((Yulab_data$Herb_pinyin_name) %in% herb), ] %>%
      dplyr::select(c(2, 4, 5)) %>%
      dplyr::distinct(.keep_all = TRUE) %>%
      tidyr::drop_na()
    colnames(result) <- c("herb", "molecule", "target")
    return(result)
  }
  
  else if (type == "Herb_en_name"){
    if (all(herb %in% (unique(Yulab_data$Herb_en_name))) == FALSE){
      herb_not_exist <- setdiff(herb, unique(Yulab_data$Herb_en_name))
      print(paste0(herb_not_exist, "doesn't/don't exist in our Yulab_datas."))
      herb <- herb[-match(herb_not_exist, herb)]
    }
    
    result <- Yulab_data[which((Yulab_data$Herb_en_name) %in% herb), ] %>%
      dplyr::select(c(2, 4, 5)) %>%
      dplyr::distinct(.keep_all = TRUE) %>%
      tidyr::drop_na()
    colnames(result) <- c("herb", "molecule", "target")
    return(result)
  }
  
  else{
    stop("parameter 'type' should be consistent of 'Herb_ch_name', 'Herb_pinyin_name', and 'Herb_en_name'")
  }
}


## 通过已知靶点找对应单体和中药函数
target_herb <- function(gene_list){
  
  library(dplyr)
  load("Yulab_data.rda")
  
  if (!all(gene_list %in% (unique(Yulab_data$target)))){
    gene_diff <- setdiff(gene_list, unique(Yulab_data$target))
    print(paste0(gene_diff, " doesn't/don't exist in the datasets."))
    gene <- gene[-match(gene_diff, gene_list)]
  }
  
  herbs_data <- Yulab_data[which(Yulab_data$target %in% gene_list), ] %>%
    dplyr::select(c(2,4,5)) %>%
    dplyr::distinct(.keep_all = TRUE) %>%
    tidyr::drop_na()
  
  colnames(herbs_data) <- c("herb", "molecule", "target")
  return(herbs_data)
}


## 对接clusterprofiler进行KEGG富集分析
enrich_target_KEGG <- function(
    x, ## 一个dataframe，包含靶点信息
    gene_list = x$target, ## 一个基因集，用于与x的靶点取交集，默认为x的所有靶点
    organism = "hsa",
    OrgDb = "org.Hs.eg.db",
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    qvalueCutoff = 0.2
){
  
  if (!"target" %in% colnames(x)){
    stop("x should contain target information")
  }
  
  gene_input <- intersect(x$target, gene_list)
  
  if(length(gene_input) == 0){
    stop("no intersection between two gene-set")
  }
  
  eg <- bitr(unique(gene_input), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = OrgDb)
  
  enrich_result <- enrichKEGG(
    gene = eg$ENTREZID,
    organism = organism,
    pvalueCutoff = pvalueCutoff,
    pAdjustMethod = pAdjustMethod,
    qvalueCutoff = qvalueCutoff
  )
  enrich_result <- setReadable(enrich_result, "org.Hs.eg.db", keyType = "ENTREZID")
  
  return(enrich_result)
  
}

## 对接clusterprofiler进行GO富集分析
enrich_target_GO <- function(
    x, ## 一个dataframe，包含靶点信息
    gene_list = x$target, ## 一个基因集，用于与x的靶点取交集，默认为x的所有靶点
    ont = "BP",
    OrgDb = "org.Hs.eg.db",
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    qvalueCutoff = 0.2,
    readable = TRUE
){
  if (!"target" %in% colnames(x)){
    stop("x should contain target information")
  }
  
  gene_input <- intersect(x$target, gene_list)
  
  if(length(gene_input) == 0){
    stop("no intersection between two gene-set")
  }
  
  eg <- bitr(unique(gene_input), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = OrgDb)
  
  enrich_result <- enrichGO(
    gene = eg$ENTREZID,
    "org.Hs.eg.db",
    ont = ont,
    pvalueCutoff = pvalueCutoff,
    pAdjustMethod = pAdjustMethod,
    qvalueCutoff = qvalueCutoff
  )
  
  return(enrich_result)
}
