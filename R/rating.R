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
#ratingArray <- read.csv("./misc/Data/ratings.csv")

searchGraph <- function(graph){
  traversals <-  list()
  for(node in graph@nodes){
    traversal <-  depthFirstSearch(graph, node@ID)
    traversals <-  append(traversals, traversal)
  }
  return(traversals)
}

# extractAllPaths <- function(graph, dfs_traversal){
#   paths <-  list()
#   current_path <-  list()
#   prev_node <-  dfs_traversal[[1]]@source_node_id
#   curr_node <-  NULL
#   origin_node <- prev_node
#   for(edge in dfs_traversal){
#     curr_node <-  edge@source_node_id
#     if(curr_node == prev_node){
#       current_path <- append(current_path, edge@source_breakpoint_id)
#       current_path <- append(current_path, edge@destination_breakpoint_id)
#       prev_node <-  edge@destination_node_id
#     }
#     else{
#       paths <- c(paths, list(current_path))
#       search <- paste0(curr_node,"_")
#       continue_index <- grep(search, rev(current_path), fixed = T)
#       #continue_index <- grep(search, current_path)
#       if (length(continue_index) > 0) {
#         index <- length(current_path) - continue_index
#         current_path <- current_path[seq_len(index)]
#
#       } else {
#         stop("This should find some match!!!")
#         current_path <- current_path
#       }
#       current_path <- append(current_path, edge@source_breakpoint_id)
#       current_path <- append(current_path, edge@destination_breakpoint_id)
#       prev_node <- edge@destination_node_id
#     }
#   }
#   return(paths)
# }

extractAllPaths <- function(graph, dfs_traversal) {
  paths <- list()
  current_path <- c()
  prev_node <- dfs_traversal[[1]]@source_node_id
  for (edge in dfs_traversal) {
    curr_node <- edge@source_node_id
    if (curr_node != prev_node) {
      paths <- c(paths, list(current_path))
      continue_index <- grep(paste0(curr_node, "_"), rev(current_path), fixed = TRUE)
      if (length(continue_index) > 0) {
        index <- length(current_path) - continue_index
        current_path <- current_path[seq_len(index)]
      } else {
        stop("This should find some match!!!")
      }
    }
    current_path <- c(current_path, edge@source_breakpoint_id, edge@destination_breakpoint_id)
    prev_node <- edge@destination_node_id
  }
  return(paths)
}


assignPathRating <- function(graph, path){
  ratings <- list()
  totalRating <- 0
  for (i in seq(1, length(path), by = 2)) {
    # iterate by two - start and end breakpoints
    startBreakpointId <- path[[i]]
    endBreakpointId <- path[[i + 1]]

    # find breakpoint objects in graph
    startBreakpoint <-  getBreakpointByID(startBreakpointId)
    endBreakpoint <-  getBreakpointByID(endBreakpointId)

    if(class(startBreakpoint) != "Breakpoint" | class(endBreakpoint) != "Breakpoint"){
      stop("One of the breakpoints not found: ",i)
    }
    # extract Encode and Region values
    startEncode <-  startBreakpoint@encode
    startRegion <-  startBreakpoint@region
    endEncode <-  endBreakpoint@encode
    endRegion <-  endBreakpoint@region

    # find the rating
    rating <- ratingArray[startEncode, endEncode]
    ratings <- c(ratings, rating)
    totalRating <- totalRating + rating
  }
  return(totalRating)
}

analyzeGeneChains <- function(source_data){
  graph <- createGraphFromDataframe(source_data)
  traversals <- searchGraph(graph)
  allPaths <- list()
  allRatings <- c()
  for(traversal in traversals){
    paths <- extractAllPaths(graph, traversal)
    allPaths <- c(allPaths, list(paths))
  }
  for(path in allPaths){
    pathRating <- assignPathRating(graph, path)
    allRatings <- c(allRatings, pathRating)
  }
  maxRatingPath <- allPaths[[which.max(allRatings)]]
  print(maxRatingPath)
}
