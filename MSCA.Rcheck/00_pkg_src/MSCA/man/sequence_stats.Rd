% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/seq_stats.R
\name{sequence_stats}
\alias{sequence_stats}
\title{Compute sequence statistics}
\usage{
sequence_stats(
  seq_list,
  min_seq_freq = 0.01,
  min_conditional_prob = 0,
  min_relative_risk = 0
)
}
\arguments{
\item{seq_list}{A list of data frames containing sequences, typically the output of \code{\link{get_cluster_sequences}}.}

\item{min_seq_freq}{Numeric threshold (default = 0.01). Filters out sequences with relative frequency below this value.}

\item{min_conditional_prob}{Numeric threshold (default = 0). Applies only for pairwise sequences (\code{k = 2}).}

\item{min_relative_risk}{Numeric threshold (default = 0). Applies only for pairwise sequences (\code{k = 2}).}
}
\value{
A list of data frames, each containing the sequence statistics for one cluster.
}
\description{
Computes descriptive statistics for sequences, including sequence frequency for any sequence length,
and conditional probability and relative risk for sequences of length 2 (pairwise transitions).
}
\details{
For \code{k = 2}, the function computes:
\itemize{
\item \strong{seq_freq:} Proportion of all sequences that match the pair
\item \strong{conditional_prob:} P(to | from)
\item \strong{relative_risk:} conditional probability divided by the marginal probability of \code{to}
}

For \code{k > 2}, only \code{seq_freq} is computed.
}
\seealso{
\code{get_cluster_sequences}
}
