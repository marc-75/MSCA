#' Compute sequence statistics
#'
#' Computes descriptive statistics for sequences, including sequence frequency for any sequence length,
#' and conditional probability and relative risk for sequences of length 2 (pairwise transitions).
#'
#' @param seq_list A list of data frames containing sequences, typically the output of \code{\link{get_cluster_sequences}}.
#' @param min_seq_freq Numeric threshold (default = 0.01). Filters out sequences with relative frequency below this value.
#' @param min_conditional_prob Numeric threshold (default = 0). Applies only for pairwise sequences (\code{k = 2}).
#' @param min_relative_risk Numeric threshold (default = 0). Applies only for pairwise sequences (\code{k = 2}).
#'
#' @details For \code{k = 2}, the function computes:
#' \itemize{
#'   \item \strong{seq_freq:} Proportion of all sequences that match the pair
#'   \item \strong{conditional_prob:} P(to | from)
#'   \item \strong{relative_risk:} conditional probability divided by the marginal probability of \code{to}
#' }
#'
#' For \code{k > 2}, only \code{seq_freq} is computed.
#'
#' @seealso \code{get_cluster_sequences}
#' @return A list of data frames, each containing the sequence statistics for one cluster.
#' @importFrom dplyr count mutate group_by ungroup arrange rename summarise across left_join filter select
#' @importFrom dplyr everything
#' @importFrom dplyr desc

#' @export
sequence_stats <- function(seq_list,
                           min_seq_freq = 0.01,
                           min_conditional_prob = 0,
                           min_relative_risk = 0) {
  lapply(seq_list, function(df) {
    if (is.null(df) || !is.data.frame(df) || nrow(df) == 0) return(NULL)

    k <- ncol(df)
    colnames(df) <- paste0("V", seq_len(k))

    tab <- df %>%
      count(across(everything()), name = "seq_count")

    total <- sum(tab$seq_count)
    tab <- tab %>%
      mutate(seq_freq = seq_count / total)

    if (k == 2) {
      tab <- tab %>%
        rename(from = V1, to = V2) %>%
        group_by(from) %>%
        mutate(from_count = sum(seq_count),
               conditional_prob = seq_count / from_count) %>%
        ungroup()

      marginal <- tab %>%
        group_by(to) %>%
        summarise(to_prob = sum(seq_count) / total, .groups = "drop")

      tab <- tab %>%
        left_join(marginal, by = "to") %>%
        mutate(relative_risk = conditional_prob / to_prob) %>%
        select(from, to, seq_count, seq_freq, conditional_prob, relative_risk) %>%
        filter(seq_freq >= min_seq_freq,
               conditional_prob >= min_conditional_prob,
               relative_risk >= min_relative_risk)
    } else {
      tab <- tab %>%
        filter(seq_freq >= min_seq_freq)
    }

    tab %>%
      arrange(desc(seq_freq))
  })
}
