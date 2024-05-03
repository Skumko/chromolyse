values <- c(
  "intergenic",
  "intronic",
  "ncRNA_intronic",
  "ncRNA_exonic",
  "UTR5",
  "exonic",
  "downstream",
  "UTR3",
  "upstream",
  "upstream;downstream",
  "splicing",
  "ncRNA_splicing",
  "UTR5;UTR3",
  "exonic;splicing",
  "ncRNA_exonic;splicing",
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

ratingMatrix <- as.matrix(read.csv(system.file("data", "ratings_no_header.csv",package = "chromolyse"), header = F))

rownames(ratingMatrix) <- values
colnames(ratingMatrix) <- values

#' Assign rating to all Edges in a path list.
#'
#' @param graph a Graph object from which the path was created.
#' @param path the path to rank - a list of `Breakpoint_id`s.
#'
#' @return a total ranking of the path as one integer.
#' @export
#'
assignPathRating <- function(graph, path){
  ratings <- list()
  totalRating <- 0
  if(length(path) %% 2 != 0){
    stop("Incorrect path! Path should contain even number of breakpoint IDs.")
  }
  for (i in seq(1, length(path), by = 2)) {
    # iterate by two - start and end breakpoints
    startBreakpointId <- path[[i]]
    endBreakpointId <- path[[i + 1]]

    # find breakpoint objects in graph
    startBreakpoint <-  getBreakpointByID(graph, startBreakpointId)
    endBreakpoint <-  getBreakpointByID(graph, endBreakpointId)

    if(class(startBreakpoint) != "Breakpoint" | class(endBreakpoint) != "Breakpoint"){
      stop("One of the breakpoints not found: ",i)
    }
    # extract Encode and Region values
    startEncode <-  startBreakpoint@encode
    startRegion <-  startBreakpoint@region
    endEncode <-  endBreakpoint@encode
    endRegion <-  endBreakpoint@region

    # find the rating
    rating <- ratingMatrix[startEncode, endEncode] + ratingMatrix[startEncode, endRegion]
    rating <- rating + ratingMatrix[startRegion, endRegion] + ratingMatrix[startRegion, endEncode]
    rating <- rating / 4

    ratings <- c(ratings, rating)
    totalRating <- totalRating + rating
  }
  return(list(totalRating = totalRating, ratings = ratings))
  #return(totalRating)
}


analyseGeneChains <- function(sourceData, minMeanRating = 5, minSupportedReads = 3, minQuality = 70, clusteringAlgorithm = "ckmeans", maxK = 20, minNumTranslocations = 7, minNumChromosomes = 3, pathDepthLimit = 20){

  # clean the dataset
  cleanSourceData <- cleanSVDataset(sourceData, minSupportedReads = minSupportedReads, minQuality = minQuality)

  if(nrow(cleanSourceData)<2){
    stop("Filtering criteria too strict! Not enough samples for clustering.")

  }
  # cluster source_data
  clusteredDataset <- clusterData(cleanSourceData, algorithm = clusteringAlgorithm, maxK = maxK)

  # create a graph object
  graph <- createGraphFromDataframe(clusteredDataset)
  # run DFS for all possible starting nodes
  traversals <- searchGraph(graph, depthLimit = pathDepthLimit)
  allPaths <- list()
  allRatings <- c()
  for(traversal in traversals){
    # extract individual paths
    if(length(traversal)){
      pathList <- extractAllPaths(graph, traversal)
      allPaths <- c(allPaths, pathList)
    }
  }
  allPaths <- allPaths[unlist(lapply(allPaths, function(x) length(x) >= minNumTranslocations*2))]
  # rank paths
  for(path in allPaths){
    pathRating <- assignPathRating(graph, path)
    allRatings <- c(allRatings, list(pathRating))
  }

  allPathsRated <- mapply(function(path, ratings) {
    list(
      path = path,
      numChromosomes = getNumChromosomes(path),
      totalRating = ratings$totalRating,
      meanRating = mean(unlist(ratings$ratings)),
      maxRating = max(unlist(ratings$ratings)),
      ratings = ratings$ratings
    )
  }, allPaths, allRatings, SIMPLIFY = FALSE)

  # filter out paths which have less chromosomes involved than defined by minNumChromosome parameter
  allPathsRated <- allPathsRated[unlist(lapply(allPathsRated, function(x) x$numChromosomes >= minNumChromosomes & x$meanRating >= minMeanRating))]

  #order by totalRating descending
  ordered_indices <- order(sapply(allPathsRated, function(x) x$meanRating), decreasing = TRUE)
  allPathsRated <- allPathsRated[ordered_indices]

  mutations <- NULL
  if(length(allPathsRated)){
    for(i in 1:length(allPathsRated)){
      pathData <- createDataframeFromEvent(result$graph, allPathsRated[[i]])
      pathData$pathIndex <- i
      mutations <- rbind(mutations, pathData)
    }
    mutations <- mutations |> dplyr::filter(rating>= minMeanRating) |>
      dplyr::select(rating, pathIndex, gene.x, gene.y, encode.x, encode.y,
                    region.x, region.y, ID.x, ID.y, nodeID.x, position.x, nodeID.y, position.y) |>
      dplyr::arrange(desc(rating))

  }

  # Extract paths meeting the minimum translocations criterion
  #allPathsRated <- allPathsRated[unlist(lapply(allPathsRated, function(x) length(x$path) >= minNumTranslocations*2))]
  return(list(clusteredData = clusteredDataset, graph = graph, paths = allPathsRated, mutations = mutations))
}
