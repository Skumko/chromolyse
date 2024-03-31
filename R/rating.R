encodeValues <- c(
  "Intergenic",
  "Intronic",
  "ncRNA_intronic",
  "ncRNA_exonic",
  "UTR5",
  "Exonic",
  "Downstream",
  "UTR3",
  "Upstream",
  "Upstream;downstream",
  "Splicing",
  "ncRNA_splicing",
  "UTR5;UTR3",
  "Exonic;splicing",
  "ncRNA_exonic;splicing"
)

regionValues <- c(
  "Name=1_Active_Promoter",
  "Name=2_Weak_Promoter",
  "Name=3_Poised_Promoter",
  "Name=4_Strong_Enhancer",
  "Name=5_Strong_Enhancer",
  "Name=6_Weak_Enhancer",
  "Name=7_Weak_Enhancer",
  "Name=8_Insulator",
  "Name=9_Txn_Transition",
  "Name=10_Txn_Elongation",
  "Name=11_Weak_Txn",
  "Name=12_Repressed",
  "Name=13_Heterochrom/lo",
  "Name=14_Repetitive/CNV",
  "Name=15_Repetitive/CNV"
)

ratingArray <- array(NA, dim = c(length(encodeValues), length(encodeValues), length(regionValues)),
                      dimnames = list(encodeValues, encodeValues, regionValues))

# TODO assign path ratings


searchGraph <- function(graph){
  traversals <-  list()
  for(node in graph@nodes){
    traversal <-  depthFirstSearch(graph, node@ID)
    traversals <-  append(traversals, traversal)
  }
  return(traversals)
}

extractAllPaths <- function(graph, dfs_traversal){
  dfs_traversal <- path
  paths <-  list()
  current_path <-  list()
  prev_node <-  dfs_traversal[[1]]@source_node_id
  curr_node <-  NULL
  origin_node <- prev_node
  for(edge in dfs_traversal){
    curr_node <-  edge@source_node_id
    if(curr_node == prev_node){
      current_path <- append(current_path, edge@source_breakpoint_id)
      current_path <- append(current_path, edge@destination_breakpoint_id)
      prev_node <-  edge@destination_node_id
    }
    else{
      paths <- c(paths, list(current_path))
      search <- paste0(curr_node,"_")
      continue_index <- grep(search, rev(current_path), fixed = T)
      #continue_index <- grep(search, current_path)
      if (length(continue_index) > 0) {
        index <- length(current_path) - continue_index
        current_path <- current_path[seq_len(index)]

      } else {
        stop("This should find some match!!!")
        current_path <- current_path
      }
      current_path <- append(current_path, edge@source_breakpoint_id)
      current_path <- append(current_path, edge@destination_breakpoint_id)
      prev_node <- edge@destination_node_id
    }
  }
  return(paths)
}

assignPathRating <- function(graph, path){
  ratings <- list()
  for (i in 1:(length(path) - 1)) {
    # iterate by two - start and end breakpoints
    start_breakpoint_id <- path[i]
    end_breakpoint_id <- path[i + 1]

    # find breakpoint objects in graph
    start_breakpoint <-  getBreakpointByID(start_breakpoint_id)
    end_breakpoint <-  getBreakpointByID(end_breakpoint_id)

    if(class(start_breakpoint) != "Breakpoint" | class(end_breakpoint) != "Breakpoint"){
      stop("One of the breakpoints not found: ",i)
    }
    # extract Encode and Region values
    start_encode <-  start_breakpoint@encode
    start_region <-  start_breakpoint@region
    end_encode <-  end_breakpoint@encode
    end_region <-  end_breakpoint@region

    # find the rating
    rating <- ratingArray[start_encode, end_encode, start_region]
    ratings <- c(ratings, rating)
  }
  return(ratings)
}

analyzeGeneChains <- function(source_data){
  graph <- createGraphFromDataframe(source_data)
  traversals <- searchGraph(graph)
  allPaths <- list()
  for(traversal in traversals){
    paths <- extractAllPaths(graph, traversal)
  }
}


