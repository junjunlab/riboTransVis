% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/enzyme-digestion-sites.R
\name{digestion_site_plot}
\alias{digestion_site_plot}
\alias{digestion_site_plot,ribotrans-method}
\title{Visualize Nucleotide Motifs Around Digestion Sites}
\usage{
digestion_site_plot(object, ...)

\S4method{digestion_site_plot}{ribotrans}(
  object,
  transcript_fa = NULL,
  read_length = c(20, 35),
  digest_site_range = c(8, 8),
  type = c("end5", "end3"),
  show_frame = FALSE,
  method = c("prob", "bits"),
  merge_rep = FALSE,
  return_data = FALSE
)
}
\arguments{
\item{object}{A \code{ribotrans} object containing Ribo-seq summary information.}

\item{...}{Additional arguments (currently unused).}

\item{transcript_fa}{A path to a FASTA file of reference transcript sequences. This file is used to extract nucleotide sequences around the digestion site positions.}

\item{read_length}{A numeric vector of length 2 specifying the minimum and maximum read length to consider (default: \code{c(20, 35)}).}

\item{digest_site_range}{A numeric vector of length 2 specifying the number of nucleotides to extract upstream (5') and downstream (3') of the digestion site (default: \code{c(8, 8)}).}

\item{type}{A character string specifying which digestion end to use for motif extraction; either \code{"end5"} or \code{"end3"}.}

\item{show_frame}{Logical. Whether to split motif plots by ribosomal frame (0, 1, 2). Default is \code{FALSE}.}

\item{method}{A character string specifying the logo scaling method. Either \code{"prob"} for position probability matrices or \code{"bits"} for information content (bits). Passed to \code{ggseqlogo}.}

\item{merge_rep}{Logical. Whether to merge replicates based on \code{sample_group}. Default is \code{FALSE}.}

\item{return_data}{Logical. If \code{TRUE}, returns raw position frequency matrices (PFMs) instead of ggplot objects. Default is \code{FALSE}.}
}
\value{
Returns either:
\itemize{
\item A nested list of ggplot objects (one per sample, optionally separated by frame), for logo visualization.
\item A list of PFMs (if \code{return_data = TRUE}).
}
}
\description{
This method generates a sequence motif plot of digestion sites for ribosome profiling data.
The motif is derived from sequences flanking the 5' or 3' end of aligned reads in transcriptomic regions.
The function supports frame-specific motif analysis, customizable motif size, and layout adjustment.
}
\details{
This function filters input Ribo-seq reads based on read length and alignment annotation (using \code{mstart} and \code{mstop}),
calculates digestion site coordinates depending on the selected \code{type}, extracts sequences from the corresponding transcriptome,
and computes a position frequency matrix (PFM) of nucleotides using \code{Biostrings::consensusMatrix}.
The logos are plotted using \code{ggseqlogo}, optionally highlighting the digestion site and breaking down by reading frame.
}
\examples{
\dontrun{
# Load example ribotrans object and transcript FASTA file
digestion_site_plot(object = rb, transcript_fa = "trans.fa", type = "end5")

# Show frame-specific motifs in bit score
digestion_site_plot(object = rb, transcript_fa = "trans.fa", show_frame = TRUE, method = "bits")
}

}
