#' Breakpoint class
#'
#' @slot ID unique identifier. Format is `[Chromosome]_[Cluster]_[Position]`.
#' @slot node_id the node this breakpoint belongs to. See \link[chromolyse]{Node}.
#' @slot position the position of breakpoint on the chromosome.
#' @slot gene a string representing which gene the breakpoint is on.
#' @slot encode the encode value of the breakpoint.
#' @slot region the region value of the breakpoint. May also be called class in SV data.
#' @export
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
  .Object@node_id <- as.character(node_id)
  .Object@position <- as.integer(position)
  .Object@gene <- as.character(gene)
  .Object@encode <- as.character(encode)
  .Object@region <- as.character(region)
  .Object
})

#' Node class
#'
#'Represents a node in the graph network.
#' @slot ID unique identifier. Format is `[chromosome]_[cluster]`.
#' @slot chromosome which chromosome the node is on.
#' @slot cluster the cluster label assigned.
#' @slot breakpoints a list of `Breakpoints` belonging to this node.
#' @slot neighbours a list of `Breakpoints` that are the destination breakpoint in translocations which contain this node's translocations as start breakpoints.
#' @export
setClass("Node",
         representation(
           ID = "character",
           chromosome = "character",
           cluster = "integer",
           breakpoints = "list",
           neighbours = "list"
         )
)
setMethod("initialize", "Node", function(.Object, chromosome, cluster) {
  .Object@ID <- paste(chromosome, cluster, sep = "_")
  .Object@chromosome <- as.character(chromosome)
  .Object@cluster <- as.integer(cluster)
  .Object@breakpoints <- list()
  .Object@neighbours <- list()
  .Object
})

setGeneric("addBreakpoint", function(object, value) standardGeneric("addBreakpoint"))
setMethod("addBreakpoint", "Node", function(object, value) {

  if(!is(value, "Breakpoint")){
    stop("Value must be of type Breakpoint!")
  }
  object@breakpoints[[length(object@breakpoints) + 1]] <- value
  object
})

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

#' Edge class
#'
#'Represents an edge in the graph network. Also represents the translocation between two `Breakpoint` objects.
#' @slot source_breakpoint_id `Breakpoint` ID of the start breakpoint in the translocation.
#' @slot source_node_id `Node` ID of the start breakpoint in the translocation.
#' @slot source_chromosome chromosome which the start breakpoint belongs to.
#' @slot source_cluster cluster which the start breakpoint has been assigned.
#' @slot source_position start breakpoint position.
#' @slot destination_breakpoint_id `Breakpoint` ID of the end breakpoint in the translocation.
#' @slot destination_node_id ID of the end breakpoint in the translocation.
#' @slot destination_chromosome chromosome which the end breakpoint belongs to.
#' @slot destination_cluster cluster which the end breakpoint has been assigned.
#' @slot destination_position end breakpoint position.
#' @slot edge_type type of translocation: [`ITX`,`CTX`].
#' @export
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


#' Check if `Edge` is in a list of edges.
#'
#' @param edge_list the list of edges.
#' @param target_edge the `Edge` to search for.
#'
#' @return boolean whether the `Edge` was found.
#' @export
containsEdge <- function(edge_list, target_edge){
  contains_edge <- FALSE
  for(edge in edge_list){
    if(edge@source_breakpoint_id == target_edge@source_breakpoint_id &
       edge@destination_breakpoint_id == target_edge@destination_breakpoint_id &
       edge@edge_type == target_edge@edge_type){
      contains_edge <- TRUE
      break
    }
  }
  return(contains_edge)
}

#' Graph class
#'
#' Represents the graph network as a list of nodes and a list of edges.
#'
#' @slot nodes list of `Node` objects.
#' @slot edges list of `Edge` objects.
#' @export
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

setGeneric("addEdge", function(object, edge) standardGeneric("addEdge"))
setMethod("addEdge", "Graph", function(object, edge) {
  # Check if the edge is of class Edge
  if (!is(edge, "Edge")) {
    stop("Argument 'edge' must be of class 'Edge'")
  }

  # Add the edge to the list of edges in the Graph
  object@edges[[length(object@edges) + 1]] <- edge
  object
})

setGeneric("hasNode", function(object, node_id) standardGeneric("hasNode"))
setMethod("hasNode", "Graph", function(object, node_id) {
  found_node <- NULL
  for (node in object@nodes) {
    if (node@ID == node_id) {
      found_node <- node
      break
    }
  }
  return(found_node)
})

setGeneric("updateNode", function(object, node_id, new_node) standardGeneric("updateNode"))
setMethod("updateNode", "Graph", function(object, node_id, new_node) {
  index <- which(sapply(graph, function(x) x@ID == node_id))
    if (length(index) > 0) {
      object[[index]] <- new_node
    } else {
      stop("Error: Node not found in the graph")
    }
})

#' Read input data and create a `Graph` object from it.
#'
#' @param dataframe the input dataframe of translocations.
#' @return a `Graph` object
#' @export
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
    if (is.null(x_existing)) {
      new_x_node = new("Node", chromosome = x_chromosome, cluster = x_cluster)
      new_x_node <- addBreakpoint(new_x_node, x_breakpoint)
      # add neighbour only to the X (source) node since we have digraph
      new_x_node <- addNeighbour(new_x_node, y_breakpoint_id)
      graph <- addNode(graph, new_x_node)
    }
    else{
      x_existing@breakpoints[[length(x_existing@breakpoints) + 1]] <- x_breakpoint
      x_existing <- addNeighbour(x_existing, y_breakpoint_id)
      for (i in seq_along(graph@nodes)) {
        if (graph@nodes[[i]]@ID == x_existing@ID) {
          graph@nodes[[i]] <- x_existing
          break
        }
      }
    }

    y_existing <- hasNode(graph, y_node_id)
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


#' Retrieve a `Breakpoint` object by its `ID`
#'
#' @param graph the graph in which the `Breakpoint` should exist.
#' @param breakpoint the `Breakpoint` ID. Format is `[Chromosome]_[Cluster]_[Position]`.
#' @return the matching `Breakpoint` object if exists else `NULL`.
#' @export
getBreakpointByID <- function(graph, breakpoint){
  if (!is.character(breakpoint)) {
    warning("Invalid input. Expected a character.")
    return(NULL)
  }

  split_id <- strsplit(breakpoint, "_")[[1]]
  node_id <- paste(split_id[1], split_id[2], sep = "_")

  node <- NULL
  for (n in graph@nodes) {
    if (n@ID == node_id) {
      node <- n
      break
    }
  }

  if (!is.null(node)) {
    for (b in node@breakpoints) {
      if (b@ID == breakpoint) {
        return(b)
        break
      }
    }
    warning("Breakpoint not found for ID ", breakpoint)
  } else {
    warning("Node not found for ID ", node_id)
  }
}

#' Perform Depth-First Search.
#'
#' @param graph `Graph` object to search.
#' @param curr_node_id `Node` ID to start from. The algorithm is recursive and will use this parameter to modify the current node.
#' @param curr_path a list of `Edge` objects representing traversal progress. Should stay NULL when calling.
#' @param depth the depth of traversal. Should stay 0 when calling.
#' @param depthLimit an upper boundary restricting the traversal depth.
#' @return a list of `Edge` objects in order of traversal. This can be used to extract individual paths using the `extractAllPaths` function.
#' @export

depthFirstSearch <- function(graph, curr_node_id, curr_path = NULL, depth = 0, depthLimit = NULL) {

  # initialize path if empty
  if (is.null(curr_path)) {
    curr_path <- list()
  }

  if (!is.null(depthLimit) && depth >= depthLimit) {
    return(curr_path)
  }

  # extract all CTX edges
  ctx_neighbour_edges <- graph@edges[sapply(graph@edges, function(x) {
    x@source_node_id == curr_node_id && x@edge_type == "CTX"
  })]

  # extract all ITX edges
  itx_neighbour_edges <- graph@edges[sapply(graph@edges, function(x) {
    x@source_node_id == curr_node_id && x@edge_type == "ITX"
  })]

  # define the recursive approach for traversing edges
  process_edges <- function(edges) {
    for (edge in edges) {
      # if the current edge has been visited, skip
      if (containsEdge(curr_path, edge)) next

      # add edge to paths
      curr_path <- c(curr_path, edge)

      # update depth
      new_depth <- depth + 1
      # recursion
      curr_path <- depthFirstSearch(graph, edge@destination_node_id, curr_path, depth = new_depth, depthLimit = depthLimit)
    }
    return(curr_path)
  }

  curr_path <- process_edges(ctx_neighbour_edges)
  curr_path <- process_edges(itx_neighbour_edges)
  return(curr_path)
}

#' Extract exact paths from a traversal list generated by Depth-First Search (see \link[chromolyse]{depthFirstSearch}).
#'
#' The function returns a list of individual paths. The list contains Breakpoint IDs, always beginning with the starting node and organized in pairs - odd indexes are starting positions of an edge, even positions are end positions of an edge.
#' The length of such list of paths should hence always be even.
#' @param graph the graph on which the search has been performed.
#' @param dfs_traversal a sequence (list) of `Edge` objects returned by calling the `depthFirstSearch`.
#'
#' @return a list of all paths.
#' @export
extractAllPaths <- function(graph, dfs_traversal) {
  paths <- list()
  current_path <- c()
  prev_node <- dfs_traversal[[1]]@source_node_id
  for (edge in dfs_traversal) {
    curr_node <- edge@source_node_id
    if (curr_node != prev_node) {
      if(length(current_path) %% 2 == 0 ){
        # add only valid path - has even number of breakpoints
        paths <- c(paths, list(current_path))
      }
      # get the last occurrence of the node we want to continue from (in a reversed list it is first occurrence)
      continue_match <- grep(paste0(curr_node, "_"), rev(current_path), fixed = TRUE)
      if (length(continue_match) > 0) {
        continue_index <- continue_match[[1]]
        index <- length(current_path) - continue_index
        current_path <- current_path[seq_len(index)]
      } else {
        stop("Graph traversal is corrupt! Did not find node.")
      }
    }
    current_path <- c(current_path, edge@source_breakpoint_id, edge@destination_breakpoint_id)
    prev_node <- edge@destination_node_id
  }
  if(!length(paths)){
    return(list(current_path))
  }
  else{
    return(paths)
  }
}

#' Perform Depth-first search from all nodes in a graph
#'
#' @param graph the graph to search.
#'
#' @return a list of "traversals", which is a list of `Edge` objects in order they were traversed, for all nodes. See \link[chromolyse]{depthFirstSearch}.
#' @export
searchGraph <- function(graph, depthLimit = 20){
  traversals <-  list()
  for(node in graph@nodes){
    traversal <-  depthFirstSearch(graph, node@ID, depthLimit = depthLimit)
    traversals <-  append(traversals, list(traversal))
  }
  return(traversals)
}

#' Extract information from identified paths to a data table.
#'
#' @param graph the graph used for event identification
#' @param event the event object = a list generated by the `\link[chromolyse]{analyseGeneChains}` function.
#'
#' @return a data table containing all translocations present in the event.
#' @export
createDataframeFromEvent <- function(graph, event){

  columnNames <- c("ID.x","nodeID.x","position.x","gene.x","encode.x","region.x",
                   "ID.y","nodeID.y","position.y","gene.y","encode.y","region.y",
                   "type")
  df <- data.frame(matrix(ncol = length(columnNames), nrow = 0))
  colnames(df) <- columnNames

  path <- event$path

  slotNames <- slotNames(getBreakpointByID(graph, path[[1]]))

  for (i in seq(1, length(path), by = 2)) {
    originBreakpoint <- getBreakpointByID(graph, path[[i]])
    destinationBreakpoint <- getBreakpointByID(graph, path[[i+1]])
    originValues <- lapply(slotNames, function(name){
      slot(originBreakpoint, name = name)
    })
    destinationValues <- lapply(slotNames, function(name){
      slot(destinationBreakpoint, name = name)
    })
    row <- c(originValues, destinationValues)
    if(unlist(strsplit(originValues[[1]], "_"))[1] == unlist(strsplit(destinationValues[[1]], "_"))[1]){
      row <- c(row, "ITX")
    }
    else{
      row <- c(row, "CTX")
    }
    names(row) <- columnNames
    df <- rbind(df, row)

  }
  df$rating <- unlist(event$ratings)
  return(df)

}

#' Transform a list of `Breakpoint` object IDs into a dataframe.
#'
#' @param graph a `Graph` object. Used to extract `Breakpoint` objects by their IDs.
#' @param path a list of `Breakpoint` object IDs.
#'
#' @return a data table with all translocations and their respective attributes.
#' @export
createDataframeFromPath <- function(graph, path){

  columnNames <- c("ID.x","nodeID.x","position.x","gene.x","encode.x","region.x",
                   "ID.y","nodeID.y","position.y","gene.y","encode.y","region.y",
                   "type")
  df <- data.frame(matrix(ncol = length(columnNames), nrow = 0))
  colnames(df) <- columnNames

  slotNames <- slotNames(getBreakpointByID(graph, path[[1]]))

  for (i in seq(1, length(path), by = 2)) {
    originBreakpoint <- getBreakpointByID(graph, path[[i]])
    destinationBreakpoint <- getBreakpointByID(graph, path[[i+1]])
    originValues <- lapply(slotNames, function(name){
      slot(originBreakpoint, name = name)
    })
    destinationValues <- lapply(slotNames, function(name){
      slot(destinationBreakpoint, name = name)
    })
    row <- c(originValues, destinationValues)
    if(unlist(strsplit(originValues[[1]], "_"))[1] == unlist(strsplit(destinationValues[[1]], "_"))[1]){
      row <- c(row, "ITX")
    }
    else{
      row <- c(row, "CTX")
    }
    names(row) <- columnNames
    df <- rbind(df, row)

  }
  return(df)

}

#' Extract number of chromosomes from a list of `Breakpoint` IDs representing a path.
#'
#' @param path a list of `Breakpoint` IDs. Likely a result of `extractAllPaths` function.
#'
#' @return an integer representing the number of chromosomes the `path` involves.
#' @export
getNumChromosomes <- function(path) {
  # Function to extract chromosome number from breakpoint identification string
  extract_chromosome <- function(string) {
    chromosome <- strsplit(string, "_")[[1]][1]
    return(chromosome)
  }

  # Extract chromosome numbers from each breakpoint identification string
  chromosomes <- sapply(path, extract_chromosome)

  # Count the number of unique chromosome numbers
  num_unique_chromosomes <- length(unique(chromosomes))

  return(num_unique_chromosomes)
}


#' Transform `Graph` nodes into an adjancency matrix.
#'
#' The function iterates over the list of `Edge` objects in `graph@edges` and accumulates the total rating of all translocations to a combination of start and destination `Node` IDs.
#' @param graph the `Graph` object from which the matrix will be created.
#'
#' @return the adjacency matrix. Each cell has the sum of all translocation ratings.
#' @export
getAdjacencyMatrix <- function(graph){
  nodeNames <- lapply(graph@nodes, function(x) x@ID)
  adj_matrix <- matrix(0, nrow = length(graph@nodes), ncol = length(graph@nodes), dimnames = list(nodeNames, nodeNames))
  for (edge in graph@edges) {
    rating <- chromolyse::rateEdge(graph, edge@source_breakpoint_id, edge@destination_breakpoint_id)
    adj_matrix[edge@source_node_id, edge@destination_node_id] <- adj_matrix[edge@source_node_id, edge@destination_node_id] + rating
  }
  return(adj_matrix)
}
