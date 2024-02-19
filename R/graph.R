# Define Node class using S3
#' Title
#'
#' @param id
#'
#' @return
#' @export
#'
#' @examples
Node <- function(id) {
  neighbors <- list()
  class(node) <- "Node"
  node <- list(id = id, neighbors = neighbors)
  return(node)
}

# Define Edge class using S3
#' Title
#'
#' @param from
#' @param to
#' @param info
#'
#' @return
#' @export
#'
#' @examples
Edge <- function(from, to, info) {
  class(edge) <- "Edge"
  edge <- list(from = from, to = to, info = info)
  return(edge)
}

# Define Graph class using S3
#' Title
#'
#' @param directed
#'
#' @return
#' @export
#'
#' @examples
Graph <- function(directed = FALSE) {
  nodes <- list()
  class(graph) <- "Graph"
  graph <- list(nodes = nodes, edges = list(), directed = directed)
  return(graph)
}

# Add method to add edge for Node class
#' Title
#'
#' @param from
#' @param to
#' @param info
#' @param graph
#'
#' @return
#' @export
#'
#' @examples
add_edge <- function(from, to, info, graph) {
  from_id <- from$id
  to_id <- to$id
  if (!(to_id %in% names(from$neighbors))) {
    from$neighbors[[to_id]] <- TRUE
    edge <- Edge(from_id, to_id, info)
    graph$edges[[length(graph$edges) + 1]] <- edge
  } else {
    warning("Edge already exists.")
  }
}

# Add method to add node for Graph class
#' Title
#'
#' @param graph
#' @param id
#'
#' @return
#' @export
#'
#' @examples
add_node <- function(graph, id) {
  if (!(id %in% sapply(graph$nodes, function(node) node$id))) {
    node <- Node(id)
    graph$nodes[[length(graph$nodes) + 1]] <- node
  } else {
    warning("Node already exists.")
  }
}

# Function to create graph from data frame
#' Title
#'
#' @param data
#'
#' @return
#' @export
#'
#' @examples
create_graph_from_data <- function(data) {
  graph <- Graph(directed = TRUE)

  # Step 1: Sort data by Chr.x and Cluster.x columns
  data <- data[order(data$Chr.x, data$Cluster.x), ]

  # Step 2 and 3: Iterate over rows and add nodes and edges
  for (i in 1:nrow(data)) {
    # Extract information from the data frame
    chr_x <- data$Chr.x[i]
    cluster_x <- data$Cluster.x[i]
    chr_y <- data$Chr.y[i]
    cluster_y <- data$Cluster.y[i]
    position_x <- data$Position.x[i]
    position_y <- data$Position.y[i]
    encode_x <- data$Encode.x[i]
    encode_y <- data$Encode.y[i]
    gene_x <- data$Gene.x[i]
    gene_y <- data$Gene.y[i]

    # Create source and target IDs
    source_id <- paste(chr_x, cluster_x, sep = "_")
    target_id <- paste(chr_y, cluster_y, sep = "_")

    # Add nodes if they don't exist
    add_node(graph, source_id)
    add_node(graph, target_id)

    # Create edge information
    edge_info <- list(
      start_position = position_x,
      end_position = position_y,
      start_encode = encode_x,
      end_encode = encode_y,
      start_gene = gene_x,
      end_gene = gene_y
    )

    # Add edge between source and target nodes
    add_edge(graph$nodes[[source_id]], graph$nodes[[target_id]], edge_info, graph)
  }

  return(graph)
}

# Example usage
# Assuming 'data' is your data frame
# data <- read.csv("your_data.csv")

# Create graph from data
# graph <- create_graph_from_data(data)
