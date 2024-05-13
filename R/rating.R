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


#' Assing rating to an edge.
#'
#' Function assigns rating to an arbitrary edge in a graph network generated from structural variation data.
#' The edge represents a translocation, defined by origin and target breakpoints and their properties.
#'
#' @param graph The graph network object `Graph` which should contain the edge.
#' @param sourceBreakpointId The string identificator of translocation origin breakpoint.
#' @param destinationBreakpointId The string identificator of translocation target breakpoint.
#'
#' @return The rating of the edge as a floating point value.
#' @export
rateEdge <- function(graph, sourceBreakpointId, destinationBreakpointId){
  sourceBreakpoint <- chromolyse::getBreakpointByID(graph,sourceBreakpointId)
  destinationBreakpoint <- chromolyse::getBreakpointByID(graph,destinationBreakpointId)

  if(class(sourceBreakpoint) != "Breakpoint" | class(destinationBreakpoint) != "Breakpoint"){
    stop("One of the breakpoints not found: ",i)
  }
  # extract Encode and Region values
  startEncode <-  sourceBreakpoint@encode
  startRegion <-  sourceBreakpoint@region
  endEncode <-  destinationBreakpoint@encode
  endRegion <-  destinationBreakpoint@region

  # find the rating
  rating <- ratingMatrix[startEncode, endEncode] +
    ratingMatrix[startEncode, endRegion] +
    ratingMatrix[startRegion, endRegion] +
    ratingMatrix[startRegion, endEncode]
  rating <- rating / 4
  return(rating)
}

#' Assign rating to a path (list of graph breakpoints).
#'
#' @param graph A `Graph` object from which the path was created.
#' @param path The path to rank - a list of `Breakpoint_ID`s.
#'
#' @return The total ranking of the path as as floating point number.
#' @export
assignPathRating <- function(graph, path){
  ratings <- list()
  totalRating <- 0
  if(length(path) %% 2 != 0){
    stop("Incorrect path! Path should contain even number of breakpoint IDs.")
  }
  for (i in seq(1, length(path), by = 2)) {
    # iterate by two - start and end breakpoints
    sourceBreakpointId <- path[[i]]
    destinationBreakpointId <- path[[i + 1]]

    rating <- rateEdge(graph, sourceBreakpointId, destinationBreakpointId)
    ratings <- c(ratings, rating)
    totalRating <- totalRating + rating
  }
  return(list(totalRating = totalRating, ratings = ratings))
}


#' Identify structural rearrangement events from the source data.
#'
#' The function takes in the source structural variation data, performs necessary preprocessing,
#' and using graph network search identifies potential massive structural rearrangement events
#' based on the defined criteria.
#'
#' @param sourceData The structural variation dataset.
#' @param cnvData Optional CNV data used in visualisation.
#' @param minMeanRating The threshold for the minimum mean rating of identified paths. For example, setting `minMeanRating=4`,
#' the function only returns identified paths which have mean rating 4 and higher.
#' @param minSupportedReads The threshold for the minimum number of supporting reads for a translocation. Used in preprocessing, see `cleanSVDataset` function.
#' @param minQuality The threshold for the minimum quality of a translocation. Used in preprocessing, see also `cleanSVDataset` function.
#' @param clusteringAlgorithm The clustering algorithm used for clustering translocation breakpoints. One of `["ckmeans", "gmm", "affinity"]` values accepted. See `clusterData` function.
#' @param maxK The upper threshold for selecting the optimal number of clusters during breakpoint clustering. See `clusterData` function.
#' @param minNumTranslocations The threshold for the minimum number of translocations in a path. For example, setting `minNumTranslocations=10`,
#'  the function only returns identified paths which contain 10 or more translocations.
#' @param minNumChromosomes The threshold for the minimum number of chromosomes involved in a path. For example, setting `minNumChromosomes=5`,
#'  the function only returns identified paths which contain translocation between 5 or more chromosomes.
#' @param pathDepthLimit The upper boundary for the depth of a path, see `depthFirstSearch` function.
#'
#' @return A list containing the results of analysis:
#' \itemize{
#' \item `clusteredData`, the source dataset annotated with assigned cluster labels,
#' \item `graph`, a `Graph` type object generated from the source data,
#' \item `events`, a list of identified chains/paths which represent potential events of massive structural rearrangement. Each event contains:
#' \itemize{
#' \item `path`, ,
#' \item `totalRating`, ,
#' \item `meanRating`, ,
#' \item `maxRating`, ,
#' \item `numChromosomes`, ,
#' },
#' \item `mutations`, a concatenated dataset of translocations that are part of individual identified paths,
#' \item `allPaths`, a concatenated list of all paths, represented by origin and target breakpoint IDs. Used for visualisation.
#' }
#' @export
analyseGeneChains <- function(sourceData, cnvData = NULL, minMeanRating = 4.5, minSupportedReads = 3, minQuality = 70, clusteringAlgorithm = "ckmeans", maxK = 20, minNumTranslocations = 7, minNumChromosomes = 3, pathDepthLimit = NULL){

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
      pathData <- createDataframeFromEvent(graph, allPathsRated[[i]])
      pathData$pathIndex <- i
      mutations <- rbind(mutations, pathData)
    }
    mutations <- mutations |> dplyr::filter(rating >= minMeanRating) |>
      dplyr::select(rating, pathIndex, type,gene.x, gene.y, encode.x, encode.y,
                    region.x, region.y, ID.x, ID.y, nodeID.x, position.x, nodeID.y, position.y) |>
      dplyr::arrange(desc(rating))

  }

  pathVisualisation <- NULL

  if(length(allPathsRated)>0){
    pathVisualisation <- unlist(lapply(allPathsRated, function(x) x$path))
    print("Visualising all identified events...")
    if(!is.null(cnvData)){
      if(length(names(cnvData)) == 4 && all(names(cnvData) == c("chromosome", "start_position","end_position","event"))){
        visualiseCircos(sv_data = clusteredDataset, cnv_data = cnvData, path = pathVisualisation, sv_focus = "ctx")
      }
      else{
        warning("CNV data not in correct format! Skipping CNV track...")
        visualiseCircos(sv_data = clusteredDataset, path = pathVisualisation, sv_focus = "ctx")
      }
    }
    else{
      visualiseCircos(sv_data = clusteredDataset, path = pathVisualisation, sv_focus = "ctx")
    }
  }
  else{
    print("No events were identified.")
  }

  # Extract paths meeting the minimum translocations criterion
  #allPathsRated <- allPathsRated[unlist(lapply(allPathsRated, function(x) length(x$path) >= minNumTranslocations*2))]
  return(list(clusteredData = clusteredDataset, graph = graph, events = allPathsRated, mutations = mutations, allPaths = pathVisualisation))
}
