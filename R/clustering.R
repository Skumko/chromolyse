#' Perform clustering using Ckmeans 1D algorithm
#'
#'
#' @param positions explanation
#' @param k_upper explanation
#' @keywords internal
#' @return return value description
.clusterCkmeans <- function(positions, k_upper = 20){
  # here, the upper bound is set using the parameter k_upper, unless the number of values in positions is less than k_upper
  upper_bound <- min(length(unique(positions)), k_upper)
  labels <- Ckmeans.1d.dp::Ckmeans.1d.dp(positions, k=c(1:upper_bound))$cluster
  return(labels)
}

#' Perform clustering using Expectation-Maximisation algorithm
#'
#' During the iteration of the EM algorithm, the log-likelihood of the data is calculated at each step.
#' The algorithm continues iterating until the change in log-likelihood between consecutive iterations is smaller than the tolerance (`tol`).
#' If the absolute difference in log-likelihood between two consecutive iterations falls below the tolerance value, the algorithm considers that it has converged and stops iterating.
#' @param positions explanation
#' @param k_upper explanation
#' @keywords internal
#' @return return value description
.emCluster_k <- function(positions, k = 20, max_iterations = 100, tol = 1e-6){

    # Initialize parameters randomly
    means <- sample(positions, k)
    weights <- rep(1 / k, k)
    variances <- rep(var(positions), k)

    # Initialize iteration counter and log-likelihood
    iteration <- 0
    prev_log_likelihood <- -Inf

    # EM algorithm loop
    while (iteration < max_iterations) {
      # Expectation step
      responsibilities <- matrix(0, nrow = length(positions), ncol = k)
      for (i in 1:k) {
        responsibilities[, i] <- dnorm(data, mean = means[i], sd = sqrt(variances[i])) * weights[i]
      }
      responsibilities <- responsibilities / rowSums(responsibilities)

      # Maximization step
      Nk <- colSums(responsibilities)
      weights <- Nk / length(positions)
      means <- colSums(positions * responsibilities) / Nk
      variances <- colSums((positions - means)^2 * responsibilities) / Nk

      # Calculate log-likelihood
      log_likelihood <- sum(log(rowSums(responsibilities * outer(weights, dnorm(data, mean = means, sd = sqrt(variances)), "*"))))

      # Check convergence
      if (abs(log_likelihood - prev_log_likelihood) < tol) {
        break
      }

      # Update iteration counter and log-likelihood
      iteration <- iteration + 1
      prev_log_likelihood <- log_likelihood
    }

    # Return final parameters and cluster assignments
    return(list(weights = weights, means = means, variances = variances, responsibilities = responsibilities))

}

#' Cluster SV data
#'
#' Function performs clustering of SV data chromosome-wise using selected algorithm.
#' @param dataset The input data.
#' @param algorithm Algorithm used to cluster position values from the dataset.
#' @param grouping_column Name of the column which contains chromosome assignments. This column will be used for grouping.
#' @param start_position_column Name of the column which contains starting positions of translocations on the chromosome. These will be used for clustering.
#' @param end_position_column Name of the column which contains end positions of translocations on the chromosome. These will be used for clustering.
#' @return original dataset with cluster labels in Cluster.x column for clustered Chr.x positions and cluster labels in Cluster.y for Chr.y positions
#' @export
clusterData <- function(dataset, algorithm='ckmeans', max_k = 20){

  combined_data <- rbind(
    dataset |> dplyr::select(Chr.x, Position.x) |> dplyr::rename(Chromosome = Chr.x, Position = Position.x) |> dplyr::mutate(Type = "X"),
    dataset |> dplyr::select(Chr.y, Position.y) |> dplyr::rename(Chromosome = Chr.y, Position = Position.y) |> dplyr::mutate(Type = "Y")
  )

  switch (algorithm,
          'ckmeans' = {

            result <- combined_data |>
              dplyr::group_by(Chromosome) |>
              dplyr::mutate(Cluster = .clusterCkmeans(Position, k_upper = max_k))

          },
          'kmeans++' = {

          },
          'expectation-max' = {

            result <- combined_data |>
              dplyr::group_by(Chromosome) |>
              dplyr::mutate(Cluster = .clusterCkmeans(Position, k_upper = max_k))

          },
          'default' = {
            "Invalid algorithm selected!"
            return(NULL)
          }
  )

  merge_x <- result |> dplyr::filter(Type == 'X') |> dplyr::select(Chromosome, Position, Cluster)
  merge_y <- result |> dplyr::filter(Type == 'Y') |> dplyr::select(Chromosome, Position, Cluster)

  merged_df <- merge(df, merge_y, by.x = c("Chr.y","Position.y"), by.y = c("Chromosome","Position"), all.x = TRUE) |>
    dplyr::rename(Cluster.y = Cluster) |>
    merge(merge_x, by.x = c("Chr.x","Position.x"), by.y = c("Chromosome","Position"), all.x = TRUE) |>
    dplyr::rename(Cluster.x = Cluster) |>
    dplyr::select(ID, Chr.x, Cluster.x, Position.x, Chr.y, Cluster.y, Position.y, everything())
  return(merged_df)
}

