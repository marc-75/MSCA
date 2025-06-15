#' Fast CLARA-like clustering using Jaccard dissimilarity
#'
#' Implements a CLARA (Clustering Large Applications) strategy using Jaccard dissimilarity
#' computed on individual patients state matrices. The algorithm repeatedly samples subsets of the data,
#' performs PAM or CLARANS partition on each subset, and selects the medoids that minimise the total dissimilarity across the full dataset. Final assignments are made by mapping all data points to the nearest selected medoid.
#' Cost can also be computed on a fraction of the dataset to reduce computational cost (see details).
#'
#' @details This function adapts the original CLARA method described by Kaufman and Rousseeuw (1990)
#' in "Finding Groups in Data: An Introduction to Cluster Analysis". In addition, it allows to replace the #' PAM partitioning method by the CLARANS procedure which is less costly.
#'
#' @references
#' Kaufman, L. & Rousseeuw, P. J. (1990). *Finding Groups in Data: An Introduction to Cluster Analysis*. Wiley.
#'
#' @param data A state matrix of censored time-to-event indicators as computed by the \code{make_state_matrix} function.
#' @param k Number of returned clusters.
#' @param samples Number of random samples drawn from the analysed population.
#' @param samplesize Number of patients per sample (default: min(50 + 5k, ncol(data))).
#' @param part_method Partition method used in the partition phase `pam` (default) or `clarans`
#' @param maxneighbor If part_method = `clarans`, set the number of medoid swap (default 0,075 * samplesize).
#' @param numlocal If part_method = `clarans`, set the number of medoid swap phases (default 3).
#' @param seed Random seed for reproducibility (default: 123).
#' @param cost_comp_ratio Proportion of data sampled to compute the clustering cost (default: 1).
#' @importFrom fastkmedoids fastpam fastclarans
#' @return A list with index of patients from the sample a, medoid indices, cluster assignment, and cost.
#' \describe{
#'   \item{clustering}{An integer vector of cluster assignments for each patient.}
#'   \item{medoids}{Indices of medoids associated witht the lower cost.}
#'   \item{sample}{Indices of the sampled columns used in clustering.}
#'   \item{cost}{Total cost (median of dissimilarities to assigned medoids).}
#' }
#' @note
#' To improve efficiency, the function uses either the \code{fastpam} or \code{fastclarans}
#' procedure from the \pkg{fastkmedoids} package along with an optimised Jaccard index computation.
#'
#' For testing purposes, the \code{cost_comp_ratio} parameter allows cost evaluation
#' to be performed on a subsample of the dataset at each iteration. The final clustering cost
#' is always computed on the full dataset using the medoids that achieved the lowest cost
#' on the subsampled data.
#'
#' Use this parameter with caution: setting \code{cost_comp_ratio} < 1 introduces randomness
#' into cost estimation during the iterative phase. For final analyses, it is recommended
#' to set \code{cost_comp_ratio = 1} to ensure a deterministic cost evaluation, as implemented in the
#' original CLARA method.

#' @export
fast_clara_jaccard <- function(data, k, samples = 20, samplesize = NULL,
                               seed = 123 ,
                               cost_comp_ratio = 1 ,
                               part_method = 'pam' ,
                               maxneighbor = 0.075 ,
                               numlocal = 3 ) {
### match args
  tryCatch({
    part_method <- match.arg(part_method, choices = c("pam", "clarans"))
  }, error = function(e) {
    stop("Invalid 'part_method': must be 'pam' or 'clarans'.", call. = FALSE)
  })
###
  if (cost_comp_ratio < 0.1) {
    message("cost_comp_ratio set to 0.1 for more reliable cost estimation.")
    cost_comp_ratio <- 0.1
  }
###
  if (cost_comp_ratio > 1) {
    message("cost_comp_ratio set to 1.")
    cost_comp_ratio <- 1
  }
###
    if (inherits(data, "dist"))
    data <- as.matrix( data )
  set.seed( seed )
  n <- ncol( data )
  if (is.null(samplesize)) samplesize <- min(40 + 3 * k, n)
  ###
  if (  samplesize / k <  2 ) {
    stop("Ratio samplesize / k <  2. Consider increasing  the sample sise or decreasing k.")
  }
  ###

  nc <- ncol(data)
  cn <- colnames(data)
  costs <- vector(length = samples)
  best_cost <- Inf
  best_result <- NULL
  s_dist <- list()
  s <- 1
  for (s in seq_len(samples)) {
    message("Sample ", s)

    # Step 1: Sample indices
    cns <- sample( cn , samplesize)
    sample_data <- data[, cns , drop = FALSE]

    # Step 2: Compute Jaccard distance on sample (upper triangle)
    d_sample <- jaccard_index_rcpp_upper(sample_data)
    dimnames( d_sample ) <- list( cns , cns )
    d_sample[is.na(d_sample)] <- 1
    diag(d_sample) <- 0


    # Step 3: PAM clustering on sample
    if( part_method == 'pam' ){
      part_res <- fastkmedoids::fastpam(rdist = as.dist(t(d_sample)), n = samplesize, k = k)
    } else {
    part_res <- fastkmedoids::fastclarans(as.dist(t(d_sample)) , n = samplesize, k = k ,
                                          numlocal = numlocal,
                                          maxneighbor = maxneighbor )
    }

    # Assignment
    medoids_global <- cns[ part_res@medoids ]

    # Step 4: Extract medoid data
    medoid_data <- data[ , medoids_global , drop = FALSE]

    # Step 5: Full dissimilarity: medoids × full data
    #idx <- setdiff( seq_len( nc ) ,  medoids_global_idx )
    nc2 <- round( nc * cost_comp_ratio )
    if( cost_comp_ratio != 1){
      idx2 <- sample( cn , size =  nc2 )
    } else {
      idx2 <- cn
    }

    # Here the computation in done on the full dataset (or fractionned)
    d_full <- jaccard_index_rcpp_parallel( medoid_data , data[ , idx2 ] )
    dimnames(d_full) <- list( medoids_global , idx2 )
    diss_to_medoids <- t(d_full)
    diss_to_medoids[is.na(diss_to_medoids)] <- 1

    # Step 6: Assign each column to closest medoid
    cluster_assign_idx <- apply( diss_to_medoids , 1 , which.min )

    # Step 7: Compute cost
    cost_vector <- apply( diss_to_medoids , 1 , min )
    total_cost <- median(cost_vector)
    costs[s] <- median(cost_vector)

    # Step 8: Update best result if needed
    if (total_cost < best_cost) {
      best_cost <- total_cost
      best_result <- list(
        medoids = medoids_global,
        clustering = cluster_assign_idx,
        cost = total_cost
      )
    }
  }
  # End of loop
  # If cost_comp_ratio < 1 need to compute the full cost and assignment
  # Else nothing

  if( cost_comp_ratio < 1 ){
    medoids_global <- best_result$medoids
    # Full dissimilarity for final assignment
    medoid_data <- data[, medoids_global, drop = FALSE]
    d_full <- jaccard_index_rcpp_parallel(medoid_data, data)
    dimnames(d_full) <- list( medoids_global , cn )
    diss_to_medoids <- t(d_full)
    diss_to_medoids[is.na(diss_to_medoids)] <- 1

    # Step 6: Assign each column to closest medoid
    cluster_assign_idx <- apply( diss_to_medoids , 1 , which.min )

    # Step 7: Compute cost
    cost_vector <- apply( diss_to_medoids , 1 , min )
    total_cost <- median( cost_vector )
    ###

    # Update best_result
    best_result <- list(
      medoids = medoids_global,
      clustering = cluster_assign_idx,
      cost = total_cost
    )
  }
gc()
return(best_result)
}
