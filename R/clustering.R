#' Perform clustering using Ckmeans 1D algorithm
#'
#'
#' @param positions explanation
#' @param k_upper explanation
#' @keywords internal
#' @return return value description
.cluster_ckmeans <- function(positions, k_upper = 20){
  # here, the upper bound is set using the parameter k_upper, unless the number of values in positions is less than k_upper
  upper_bound <- min(length(unique(positions)), k_upper)
  labels <- Ckmeans.1d.dp::Ckmeans.1d.dp(positions, k=c(1:upper_bound))$cluster
  return(labels)
}

#' Cluster SV data
#'
#' Function performs clustering of SV data chromosome-wise using selected algorithm.
#' @param dataset The input data.
#' @param algorithm Algorithm used to cluster position values from the dataset.
#' @param grouping_column Name of the column which contains chromosome assignments. This column will be used for grouping.
#' @param start_position_column Name of the column which contains starting positions of translocations on the chromosome. These will be used for clustering.
#' @param end_position_column Name of the column which contains end positions of translocations on the chromosome. These will be used for clustering.
#' @return original dataset with cluster labels in Cluster.x column
#' @export
cluster_data <- function(dataset, algorithm='ckmeans', grouping_column='Chr.x', start_position_column='Position.x', end_position_column='Position.y'){

  switch (algorithm,
          'ckmeans' = {
            clustered <- dataset %>%
              dplyr::group_by_at(grouping_column) %>%
              dplyr::mutate(Cluster.x = .cluster_ckmeans(!!dplyr::sym(start_position_column)), Cluster.y = .cluster_ckmeans(!!dplyr::sym(end_position_column)));
            return(clustered)
          },
          'kmeans++' = {

          },
          'default' = {
            "Invalid algorithm selected"
          }
  )
}


