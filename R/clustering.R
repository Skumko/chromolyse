#' Perform clustering using Ckmeans 1D algorithm
#'
#'
#' @param positions explanation
#' @param param_name2 explanation
#' @param param_name3 explanation
#' @return return value description
cluster_ckmeans <- function(positions, k_upper = 20){
  # here, the upper bound is set using the parameter k_upper, unless the number of values in positions is less than k_upper
  upper_bound <- min(length(unique(positions)), k_upper)
  labels <- Ckmeans.1d.dp(positions, k=c(1:upper_bound))$cluster
  return(labels)

}

#' Cluster SV data
#'
#' Function performs clustering of SV data chromosome-wise using selected algorithm.
#' @param dataset The input data.
#' @param algorithm Algorithm used to cluster position values from the dataset.
#' @param grouping_column Name of the column which contains chromosome assignments. This column will be used for grouping.
#' @param position_column Name of the column which contains positions on the chromosome. This will be used for clustering directly.
#' @return original dataset with cluster labels in Cluster.x column
#' @examples
#' temp1 <- <function_name>(50);
#' temp2 <- <function_name>( c(50, 63, 23) );
#' @export
cluster_data <- function(dataset, algorithm='ckmeans', grouping_column='Chr.x', position_column= 'Position.x'){

  switch (algorithm,
          'ckmeans' = {
            clustered <- dataset %>% group_by_at(grouping_column) %>% mutate(Cluster.x = cluster_ckmeans(!!sym(position_column)));
            return(clustered)
          },
          'kmeans++' = {

          },
          'default' = {
            "Invalid algorithm selected"
          }
  )
}


