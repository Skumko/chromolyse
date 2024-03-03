setClass("Breakpoint",
         representation(
           ID = "character",
           node_id = "character",
           position = "integer",
           gene = "character",
           encode = "character",
           region = "character"
         )
)
setMethod("initialize", "Breakpoint", function(.Object, node_id, position, gene, encode, region) {
  .Object@ID <- paste(node_id, position, sep = "_")
  .Object@node_id <- node_id
  .Object@position <- as.integer(position)
  .Object@gene <- gene
  .Object@encode <- encode
  .Object@region <- region
  .Object
})

setClass("Node",
         representation(
           ID = "character",
           chromosome = "character",
           cluster = "integer",
           breakpoints = "list",
           neighbours = "list",
         )
)
setMethod("initialize", "Node", function(.Object, chromosome, cluster) {
  .Object@ID <- paste(chromosome, cluster, sep = "_")
  .Object@chromosome <- chromosome
  .Object@cluster <- as.integer(cluster)
  .Object@breakpoints <- list()
  .Object@neighbours <- list()
  .Object
})

# Method to add a position to a Node
setGeneric("addBreakpoint", function(object, value) standardGeneric("addBreakpoint"))
setMethod("addBreakpoint", "Node", function(object, value) {

  if(!is(value, "Breakpoint")){
    stop("Value must be of type Breakpoint!")
  }
  object@breakpoints[[length(object@breakpoints) + 1]] <- value
  object
})

# Method to add a neighbour to a Node, checking for Node type
setGeneric("addNeighbour", function(object, neighbour_id) standardGeneric("addNeighbour"))
setMethod("addNeighbour", "Node", function(object, neighbour_id) {

  # Check if node ID already exists in the neighbour list
  if (neighbour_id %in% object@neighbours) {
    #warning("Node with ID ", node@ID, " already exists, skipping addition.")
    return(object)
  }
  object@neighbours[[length(object@neighbours) + 1]] <- neighbour_id
  object
})

setClass("Edge",
         representation(
           source_breakpoint_id = "character",
           source_node_id = "character",
           source_chromosome = "character",
           source_cluster = "integer",
           source_position = "integer",
           destination_breakpoint_id = "character",
           destination_node_id = "character",
           destination_chromosome = "character",
           destination_cluster = "integer",
           destination_position = "integer",
           edge_type = "character"
         )
)
setMethod("initialize", "Edge", function(.Object, source_breakpoint_id, destination_breakpoint_id, edge_type) {
  source_parts <- strsplit(source_breakpoint_id, "_")[[1]]
  .Object@source_breakpoint_id <- source_breakpoint_id
  .Object@source_node_id <- paste(source_parts[1], source_parts[2], sep = "_")
  .Object@source_chromosome <- source_parts[1]
  .Object@source_cluster <- as.integer(source_parts[2])
  .Object@source_position <- as.integer(source_parts[3])

  destination_parts <- strsplit(destination_breakpoint_id, "_")[[1]]
  .Object@destination_breakpoint_id <- destination_breakpoint_id
  .Object@destination_node_id <- paste(destination_parts[1], destination_parts[2], sep = "_")
  .Object@destination_chromosome <- destination_parts[1]
  .Object@destination_cluster <- as.integer(destination_parts[2])
  .Object@destination_position <- as.integer(destination_parts[3])

  .Object@edge_type <- edge_type
  .Object
})


setClass("Graph",
         representation(
           nodes = "list",
           edges = "list"
         )
)
setMethod("initialize", "Graph", function(.Object, nodes = list(), edges = list()) {
  .Object@nodes <- nodes
  .Object@edges <- edges
  .Object
})

# Method to add a Node object to the Graph
setGeneric("addNode", function(object, node) standardGeneric("addNode"))
setMethod("addNode", "Graph", function(object, node) {
  # Check if the node is of class Node
  if (!is(node, "Node")) {
    stop("Argument 'node' must be of class 'Node'")
  }

  # Check if node ID already exists
  existing_node_ids <- sapply(object@nodes, function(n) n@ID)
  if (node@ID %in% existing_node_ids) {
    #warning("Node with ID ", node@ID, " already exists, skipping addition.")
    return(object)
  }

  # Add the node to the list of nodes in the Graph
  object@nodes[[length(object@nodes) + 1]] <- node
  object
})

# Method to add an Edge object to the Graph
setGeneric("addEdge", function(object, edge) standardGeneric("addEdge"))
setMethod("addEdge", "Graph", function(object, edge) {
  # Check if the edge is of class Edge
  if (!is(edge, "Edge")) {
    stop("Argument 'edge' must be of class 'Edge'")
  }

  # Check if edge already exists (considering both directions)
  #existing_edges <- lapply(object@edges, function(e) c(e@source_breakpoint_id, e@destination_breakpoint_id))
  #existing_edges <- unique(unlist(existing_edges))
  #if (length(existing_edges) > 0 & edge@source_breakpoint_id %in% existing_edges & edge@destination_breakpoint_id %in% existing_edges) {
  #  warning("Edge already exists between", edge@source_breakpoint_id, "and", edge@destination_breakpoint_id, ", skipping addition.")
  #  return(object)
  #}

  # Add the edge to the list of edges in the Graph
  object@edges[[length(object@edges) + 1]] <- edge
  object
})

# Method to check if a node with specific ID exists in the graph, returning node or FALSE
setGeneric("hasNode", function(object, node_id) standardGeneric("hasNode"))
setMethod("hasNode", "Graph", function(object, node_id) {
  # Find nodes with matching ID
  matching_nodes <- sapply(object@nodes, function(n) if (n@ID == node_id) n else NULL)
  # Return the first matching node or NULL
  if(length(matching_nodes) == 0)
    return(NULL)
  possible_node <- matching_nodes[[1]]
  return(possible_node)
})

createGraphFromDataframe <- function(dataframe) {
  # Initialize empty graph, nodes, breakpoints, and edges lists
  graph <- new("Graph", nodes = list(), edges = list())

  # Iterate through rows
  for (i in 1:nrow(dataframe)) {
    # Extract data from each row
    x_chromosome <- dataframe$Chr.x[i]
    x_cluster <- dataframe$Cluster.x[i]
    x_position <- dataframe$Position.x[i]
    x_gene <- dataframe$Gene.x[i]
    x_encode <- dataframe$Encode.x[i]
    x_region <- dataframe$Class.x[i]
    y_chromosome <- dataframe$Chr.y[i]
    y_cluster <- dataframe$Cluster.y[i]
    y_position <- dataframe$Position.y[i]
    y_gene <- dataframe$Gene.y[i]
    y_encode <- dataframe$Encode.y[i]
    y_region <- dataframe$Class.y[i]
    type <- dataframe$Type[i]

    x_node_id <- paste(x_chromosome, x_cluster, sep = "_")
    y_node_id <- paste(y_chromosome, y_cluster, sep = "_")
    x_breakpoint_id <- paste(x_node_id, x_position, sep = "_")
    y_breakpoint_id <- paste(y_node_id, y_position, sep = "_")


    # Create new objects - Breakpoints and Edge
    x_breakpoint <- new("Breakpoint", node_id = x_node_id, position = x_position, gene = x_gene, encode = x_encode, region = x_region)
    y_breakpoint <- new("Breakpoint", node_id = y_node_id, position = y_position, gene = y_gene, encode = y_encode, region = y_region)
    new_edge <- new("Edge",source_breakpoint_id = x_breakpoint_id,
                destination_breakpoint_id = y_breakpoint_id,
                edge_type = type)

    #find whether Nodes already exist in graph
    x_existing <- hasNode(graph, x_node_id)
    y_existing <- hasNode(graph, y_node_id)


    if (is.null(x_existing)) {
      new_x_node = new("Node", chromosome = x_chromosome, cluster = x_cluster)
      new_x_node <- addBreakpoint(new_x_node, x_breakpoint)
      # add neighbour only to the X (source) node since we have digraph
      new_x_node <- addNeighbour(new_x_node, y_node_id)
      graph <- addNode(graph, new_x_node)
    }
    else{
      x_existing@breakpoints[[length(x_existing@breakpoints) + 1]] <- x_breakpoint
      x_existing <- addNeighbour(x_existing, y_node_id)
      for (i in seq_along(graph@nodes)) {
        if (graph@nodes[[i]]@ID == x_existing@ID) {
          graph@nodes[[i]] <- x_existing
          break
        }
      }
    }
    if (is.null(y_existing)) {
      new_y_node = new("Node", chromosome = y_chromosome, cluster = y_cluster)
      new_y_node <- addBreakpoint(new_y_node, y_breakpoint)
      graph <- addNode(graph, new_y_node)

    }
    else{
      y_existing@breakpoints[[length(y_existing@breakpoints) + 1]] <- y_breakpoint
      for (i in seq_along(graph@nodes)) {
        if (graph@nodes[[i]]@ID == y_existing@ID) {
          graph@nodes[[i]] <- y_existing
          break
        }
      }
    }

    graph <- addEdge(graph, new_edge)
  }

  return(graph)
}


## Check for the existence of node before adding
#if ("Node1" %in% names(graph$nodes)) {
#  stop("Node with the same identifier already exists.")
#} else {
#  graph$add_node(node1)
#}
#
## Check for the existence of edge before adding
#if ("Node1" %in% sapply(graph$edges, `[[`, "source") || "Node2" %in% sapply(graph$edges, `[[`, "destination")) {
#  stop("Edge with the same source or destination already exists.")
#} else {
#  graph$add_edge(edge1)
#}

#  # Method to add a node to the graph
#  graph$add_node <- function(node) {
#    if (node$identifier %in% names(graph$nodes)) {
#      stop("Node with the same identifier already exists.")
#    } else {
#      graph$nodes[[node$identifier]] <- node
#    }
#  }
#
#  # Method to add an edge to the graph
#  graph$add_edge <- function(edge) {
#    if (edge$source %in% sapply(graph$edges, `[[`, "source") || edge$destination %in% sapply(graph$edges, `[[`, "destination")) {
#      stop("Edge with the same source or destination already exists.")
#    } else {
#      graph$edges <- c(graph$edges, list(edge))
#    }
#  }
#
#  return(graph)
#}


## Define the Depth-First Search function
#dfs <- function(graph, node, visited_positions = NULL) {
#  # Check if all positions of the current node have been visited
#  all_positions_visited <- length(node$positions) == length(visited_positions)
#
#  # If all positions have been visited, mark the node as visited
#  if (all_positions_visited) {
#    node$visited <- TRUE
#    cat("Visited node:", node$identifier, "\n")
#
#    # Traverse neighboring nodes
#    for (edge in graph$edges) {
#      if (edge$source == node$identifier) {
#        destination_node <- graph$nodes[[edge$destination]]
#
#        # Check if the destination node has any unvisited positions
#        unvisited_positions <- setdiff(destination_node$positions, visited_positions)
#
#        # If there are unvisited positions, recursively visit the node with them
#        if (length(unvisited_positions) > 0) {
#          dfs(graph, destination_node, c(visited_positions, unvisited_positions))
#        }
#      }
#    }
#  }
#}
#
## Example usage:
## Assuming you have already created the graph and added nodes and edges#
#
## Choose a starting node to begin the DFS traversal
#start_node <- graph$nodes[["Node1"]]
#
## Perform DFS traversal from the starting node
#dfs(graph, start_node)

