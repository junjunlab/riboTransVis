% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gene-expression-analysis.R
\name{get_normalized_reads}
\alias{get_normalized_reads}
\alias{get_normalized_reads,ribotrans-method}
\title{Get Normalized Read Counts (TPM) and Translation Efficiency}
\usage{
get_normalized_reads(object, ...)

\S4method{get_normalized_reads}{ribotrans}(object, type = c("rna", "ribo", "te"))
}
\arguments{
\item{object}{An S4 object of class \code{ribotrans}, which stores:
• \code{object@counts}: a list containing read count matrices for RNA-seq (\code{object@counts$rna})
and Ribo-seq (\code{object@counts$rpf}).
• \code{object@features}: a data frame or tibble that holds transcript or gene length information
for "transcriptome" mapping mode.
• \code{object@mapping_type}: a character string, either \code{"genome"} or \code{"transcriptome"}.}

\item{...}{Additional arguments (currently unused).}

\item{type}{A character vector specifying which data to process. Can be any combination of
\code{c("rna","ribo","te")}.
• "rna": compute TPM for RNA-seq.
• "ribo": compute TPM for Ribo-seq.
• "te": compute translation efficiency = ribo_TPM / rna_TPM.}
}
\value{
A named list with up to three elements (depending on the \code{type} input):
• \code{$tpm.rna}: A data frame of TPM values for RNA-seq genes (genome mapping) or transcripts (transcriptome mapping).
• \code{$tpm.ribo}: A data frame of TPM values for Ribo-seq.
• \code{$te}: A data frame of translation efficiency, computed only when both RNA and Ribo data are available and "te" is included in \code{type}.
}
\description{
This function calculates TPM (Transcripts Per Million) values for RNA-seq and/or Ribo-seq data
(depending on the chosen type), then computes translation efficiency (TE) when both are requested.
It returns a list with elements tpm.rna, tpm.ribo, and te (depending on user input).
}
\details{
\enumerate{
\item If "rna" is selected in the argument \code{type}, the function extracts the counts from
\code{object@counts$rna} and calculates TPM for RNA-seq data.
\item If "ribo" is selected in \code{type}, the function extracts counts from
\code{object@counts$rpf} (ribosome-protected fragments) to calculate TPM for Ribo-seq data.
\item If "te" (translation efficiency) is selected, the function will intersect the rownames
of the RNA and Ribo TPM matrices to compute TE as tpm.ribo / tpm.rna for the matched genes.
}

The function handles two potential mapping types:
• "genome": Genomic annotation is used (e.g., from \code{ribo_info$annotation}, \code{rna_info$annotation}).
• "transcriptome": Transcript-based annotation, where \code{object@features} provides transcript lengths
in \code{translen} and gene IDs/annotations.
}
\examples{
\dontrun{
  # Suppose we have a 'ribotrans' object named 'ribo_obj'
  # that contains counts for RNA-seq (ribo_obj@counts$rna)
  # and Ribo-seq (ribo_obj@counts$rpf).

  # Calculate only RNA TPM:
  res_rna <- get_normalized_reads(ribo_obj, type = "rna")
  head(res_rna$tpm.rna)

  # Calculate RNA and Ribo TPM plus TE:
  res_te <- get_normalized_reads(ribo_obj, type = c("rna","ribo","te"))
  head(res_te$tpm.rna)
  head(res_te$tpm.ribo)
  head(res_te$te)
}

}
