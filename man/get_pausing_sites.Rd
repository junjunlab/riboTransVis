% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pausing-site-calculation.R
\name{get_pausing_sites}
\alias{get_pausing_sites}
\alias{get_pausing_sites,ribotrans-method}
\title{Identify Ribosome Pausing Sites from Ribosome Profiling Data}
\usage{
get_pausing_sites(object, ...)

\S4method{get_pausing_sites}{ribotrans}(
  object,
  sample_selected = NULL,
  merge_rep = FALSE,
  do_offset_correct = FALSE,
  position_shift = 0,
  exclude_length = c(0, 0),
  min_counts = 64,
  window = 200
)
}
\arguments{
\item{object}{A RiboSeq object containing ribosome profiling data with
summary_info and features slots}

\item{...}{Additional arguments (currently unused).}

\item{sample_selected}{Character vector. Names of samples to include in the
analysis. If NULL (default), all samples will be processed}

\item{merge_rep}{Logical. Whether to merge technical replicates by taking
the mean of normalized values within sample groups. Default is FALSE}

\item{do_offset_correct}{Logical. Whether to perform read offset correction
to account for ribosome A-site positioning. Default is FALSE}

\item{position_shift}{Numeric. The number of nucleotides to shift reads
when offset correction is applied. Only used when do_offset_correct = TRUE.
Default is 0}

\item{exclude_length}{Numeric vector of length 2. Number of nucleotides to
exclude from the start and end of coding sequences.
Default is c(0,0)}

\item{min_counts}{Numeric. Minimum total read counts required for a gene
to be included in the analysis. Genes with fewer counts are filtered out.
Default is 64}

\item{window}{Numeric. Window size (in nucleotides) used for calculating
pause scores. The algorithm uses overlapping windows to assess local
ribosome density. Default is 200}
}
\value{
A data.table containing pausing site information with the following columns:
\describe{
\item{sample}{Sample identifier}
\item{rname}{Gene/transcript identifier}
\item{relst}{Relative position from CDS start}
\item{normsm}{Normalized ribosome density (read count / average count)}
\item{sum_1}{Sum of normalized densities in window pos, pos + window}
\item{sum_2}{Sum of normalized densities in window pos + window/2, pos + 1.5*window}
\item{pause_score}{Calculated pause score using PausePred algorithm}
\item{total_sum}{Total sum of sum_1 and sum_2}
\item{min_sum}{Minimum of sum_1 and sum_2}
\item{max_sum}{Maximum of sum_1 and sum_2}
\item{balance_score}{Balance score (min_sum / total_sum) for quality control}
\item{imbalance_ratio}{Ratio of max_sum to min_sum for quality assessment}
\item{z_score}{Z-score of pause score relative to the gene}
}
}
\description{
This function identifies ribosome pausing sites by calculating pause scores
based on the PausePred algorithm. It processes ribosome profiling data to
detect positions where ribosomes pause during translation elongation.
}
\details{
The function implements the PausePred algorithm for identifying ribosome pausing
sites with the following steps:

\enumerate{
\item \strong{Data preprocessing}: Applies optional offset correction and filters
samples and genomic regions based on user specifications
\item \strong{Normalization}: Calculates normalized ribosome density by dividing
read counts by the average count per position for each gene
\item \strong{Quality filtering}: Excludes genes with low total read counts
(< min_counts) to ensure reliable statistical inference
\item \strong{Window-based scoring}: For each position, calculates sums in two
overlapping windows and computes pause scores using the formula:
\deqn{pause\_score = \frac{window}{2} \times normsm \times \frac{sum_1 + sum_2}{sum_1 \times sum_2}}
\item \strong{Quality metrics}: Adds balance and imbalance ratio metrics for
downstream quality control and filtering
}

The pause score reflects the likelihood that a given position represents a
genuine pausing site, with higher scores indicating stronger evidence for
ribosome pausing. The algorithm is designed to identify positions where
ribosome density is elevated relative to neighboring regions.

\strong{Window Strategy:}
\itemize{
\item Window 1: position, position + window - captures downstream density
\item Window 2: position + window/2, position + 1.5*window - overlapping window
\item The overlapping design helps distinguish genuine pauses from noise
}

\strong{Quality Control Metrics:}
\itemize{
\item \code{balance_score}: Ranges from 0 to 0.5, with higher values indicating
more balanced window sums (recommended threshold: > 0.1-0.2)
\item \code{imbalance_ratio}: Ratio of larger to smaller window sum
(recommended threshold: < 5-10)
}
}
\section{Filtering Recommendations}{

For downstream analysis, consider filtering pause sites based on:
\itemize{
\item Pause score threshold (typically 10-20 for initial screening)
\item Balance score > 0.15 to exclude sites with extremely unbalanced windows
\item Imbalance ratio < 5 to remove artifacts from coverage imbalances
\item Minimum total coverage to ensure statistical reliability
}
}

\examples{
\dontrun{
# Basic usage with default parameters
pause_sites <- get_pausing_sites(ribo_object)

# Analyze specific samples with offset correction
pause_sites <- get_pausing_sites(
  object = ribo_object,
  sample_selected = c("sample1", "sample2"),
  do_offset_correct = TRUE,
  position_shift = 12
)

# Merge replicates and use stricter filtering
pause_sites <- get_pausing_sites(
  object = ribo_object,
  merge_rep = TRUE,
  exclude_length = c(30, 30),  # Exclude first/last 30 nt
  min_counts = 100,            # Higher count threshold
  window = 150                 # Smaller window size
)

# Filter high-quality pause sites
high_quality_pauses <- pause_sites[
  pause_score >= 20 &
  balance_score >= 0.15 &
  imbalance_ratio <= 5
]

# Summary statistics
summary(pause_sites$pause_score)
table(pause_sites$sample)
}



}
\seealso{
\code{\link{do_offset_correction}} for read offset correction
}
