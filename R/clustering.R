#' Perform clustering using Ckmeans 1D algorithm
#'
#'
#' @param positions input data to cluster. This should represent a subset of positions from one chromosome.
#' @param kUpperBound the upper limit for the K value (target number of clusters).
#' @keywords internal
#' @return cluster labels.
.clusterCkmeans <- function(positions, kUpperBound = 20){
  # here, the upper bound is set using the parameter kUpperBound, unless the number of values in positions is less than kUpperBound
  upper_bound <- min(length(unique(positions)), kUpperBound)
  labels <- Ckmeans.1d.dp::Ckmeans.1d.dp(positions, k=c(1:upper_bound))$cluster
  return(labels)
}


#' Perform clustering using the Expectation-Maximisation algorithm using a Gaussian Mixture Model.
#'
#' GMM used from the Mclust library, see \code{mclust::\link[mclust:Mclust]{Mclust GMM}} for details.
#' @param positions input data to cluster. This should represent a subset of positions from one chromosome.
#' @param kUpperBound the upper limit for the K value (target number of clusters).
#' @keywords internal
#' @return cluster labels.
.clusterGMM <- function(positions, kUpperBound = 20){
  model <- mclust::Mclust(positions, G = 1:kUpperBound)
  if(is.null(model)){
    return(c(1))
  }
  return(model$classification)
}


#' Perform clustering using Affinity Propagation algorithm.
#'
#' The algorithm is an implementation from the apcluster package, see \code{apcluster::\link[apcluster:apcluster]{apcluster}} for details.
#'
#' @param positions input data to cluster. This should represent a subset of positions from one chromosome.
#' @return cluster labels.
#' @export
#'
.clusterAP <- function(positions){

  if(length(positions)<2){
    return(c(1))

  }
  apModel = apcluster::apcluster(apcluster::negDistMat(r=2), x=positions)
  # Map labels in positions to index numbers
  labels <- sapply(apcluster::labels(apModel), function(label) {
    which(apModel@exemplars==label)
  })
  return(labels)
}

#' Cluster SV data
#'
#' Function performs clustering of SV data chromosome-wise using selected algorithm.
#' @param dataset The input data.
#' @param algorithm Algorithm used to cluster position values from the dataset.
#' @param maxK The maximum number of clusters allowed for clustering. User definable, default set to 20.
#' @return original dataset with cluster labels in Cluster.x column for clustered Chr.x positions and cluster labels in Cluster.y for Chr.y positions
#' @export
clusterData <- function(dataset, algorithm='ckmeans', maxK = 20){

  combined_data <- rbind(
    dataset |> dplyr::select(Chr.x, Position.x) |> dplyr::rename(Chromosome = Chr.x, Position = Position.x) |> dplyr::mutate(Type = "X"),
    dataset |> dplyr::select(Chr.y, Position.y) |> dplyr::rename(Chromosome = Chr.y, Position = Position.y) |> dplyr::mutate(Type = "Y")
  )

  switch (algorithm,
          'ckmeans' = {

            result <- combined_data |>
              dplyr::group_by(Chromosome) |>
              dplyr::mutate(Cluster = .clusterCkmeans(Position, kUpperBound = maxK))

          },
          'affinity' = {
            result <- combined_data |>
              dplyr::group_by(Chromosome) |>
              dplyr::mutate(Cluster = .clusterAP(Position))

          },
          'gmm' = {

            result <- combined_data |>
              dplyr::group_by(Chromosome) |>
              dplyr::mutate(Cluster = .clusterGMM(Position, kUpperBound = maxK))

          },
          {
            stop("Error: Wrong algorithm type selected!")
          }
  )

  merge_x <- result |> dplyr::filter(Type == 'X') |> dplyr::select(Chromosome, Position, Cluster)
  merge_y <- result |> dplyr::filter(Type == 'Y') |> dplyr::select(Chromosome, Position, Cluster)
  merge_x <- merge_x |> dplyr::rename(Cluster.x = Cluster)
  merge_y <- merge_y |> dplyr::rename(Cluster.y = Cluster)

  merged_df <- dataset |>
    merge(merge_y, by.x = c("Chr.y","Position.y"), by.y = c("Chromosome","Position")) |>
    merge(merge_x, by.x = c("Chr.x","Position.x"), by.y = c("Chromosome","Position")) |>
    dplyr::select(ID, Chr.x, Cluster.x, Position.x, Chr.y, Cluster.y, Position.y, everything()) |>
    dplyr::distinct()
  return(merged_df)
}
