#' @title Visualize Pharmacology Network
#' @description Plot herb–compound–target network using ggraph with optional central herb node.
#' @param df A data.frame with columns: herb, molecule, target
#' @param herb_center Whether to center the herb (only works if 1 unique herb)
#' @param node_colors Named color vector for node types: Herb / Molecule / Target
#' @param circle_radius Numeric vector of radii for 3 layers
#' @param max.overlaps Max label overlaps
#' @return A ggplot2 ggraph object
#' @importFrom ggplot2
#' @importFrom ggraph
#‘ @importFrom igraph V
#' @importFrom dplyr select mutate

#' @export

phar_net <- function(
    df,
    herb_center = FALSE,
    node_colors = c(Herb = "red", Molecule = "#4DAF4A", Target = "#1E90FF"),
    circle_radius = c(1, 2, 3),
    max.overlaps = 100
){
  
  #library(igraph)
  #library(ggraph)
  #library(dplyr)
  #library(ggplot2)
  
  expected_cols <- c("herb", "molecule", "target")
  if (!all(expected_cols %in% colnames(df))) {
    stop(paste0("Input dataframe must contain columns: ", paste(expected_cols, collapse = ", ")))
  }
  
  edges_1 <- df %>% dplyr::select(from = herb, to = molecule)
  edges_2 <- df %>% dplyr::select(from = molecule, to = target)
  edges_all <- bind_rows(edges_1, edges_2)
  
  herbs <- unique(df$herb)
  molecules <- unique(df$molecule)
  targets <- unique(df$target)
  
  get_circle_coords <- function(n, r) {
    theta <- seq(0, 2 * pi, length.out = n + 1)[-1]
    data.frame(x = r * cos(theta), y = r * sin(theta))
  }
  
  layout_df <- bind_rows(
    data.frame(name = herbs, get_circle_coords(length(herbs), circle_radius[1]), group = "Herb"),
    data.frame(name = molecules, get_circle_coords(length(molecules), circle_radius[2]), group = "Molecule"),
    data.frame(name = targets, get_circle_coords(length(targets), circle_radius[3]), group = "Target")
  )
  
  g <- graph_from_data_frame(edges_all, vertices = layout_df$name, directed = FALSE)
  coords <- layout_df[match(V(g)$name, layout_df$name), ]
  coords <- coords %>% mutate(x = as.numeric(x), y = as.numeric(y))
  
  if (herb_center & length(herbs) == 1) {
    coords[coords$name == herbs[1], c("x", "y")] <- c(0, 0)
  }
  
  p <- ggraph(g, layout = "manual", x = coords$x, y = coords$y) +
    geom_edge_link(alpha = 0.4, color = "grey40") +
    geom_node_point(aes(color = coords$group), size = 5) +
    geom_node_text(aes(label = name), repel = TRUE, size = 3, max.overlaps = max.overlaps) +
    scale_color_manual(values = node_colors, name = NULL) +
    theme_void()
  
  return(p)
}
