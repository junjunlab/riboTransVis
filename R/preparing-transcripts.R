#' Extract Sequence from Genome Based on Transcript ID
#'
#' This function extracts the sequence of a transcript from the genome,
#' based on its specified in the GTF file.
#'
#' @param gtf_file The path to the GTF file containing transcript information.
#' @param genome_file The path to the genome file where the sequence is to be extracted
#' @param transcript_id The transcript ID for the transcript to be extracted.
#' @param new_id The new ID to assign to the extracted sequence.
#' @param type The type of sequence to extract ('CDS', 'exon', or '3/5UTR').
#' @param out_file The output file to write the extracted sequence to.
#' @param pythonPath The path of python to be used for reticulate.
#'
#' @return A FASTA format string of the extracted sequence.
#'
#'
#' @examples
#'\dontrun{
#' pyExtractSeq(gtf_file = "transcripts.gtf", genome_file = "genome.fa",
#' transcript_id = "ENST00000456328", type = "cds", out_file = "extracted_seq.fasta")
#'}
#' @export
pyExtractSeq <- function(gtf_file = NULL,
                         genome_file = NULL,
                         transcript_id = NULL,
                         new_id = NULL,
                         type = NULL,
                         out_file = NULL,
                         pythonPath = NULL){
  reticulate::py_config()
  if(reticulate::py_available() == FALSE){
    message("Please install python first!")
  }else{
    if(!is.null(pythonPath)){
      reticulate::use_python(pythonPath)
    }

    # check pyfaidx
    if (!reticulate::py_module_available("pyfaidx")) {
      cat("Installing pyfaidx ...\n")
      reticulate::py_install("pyfaidx")
    }

    # import pyfaidx
    tryCatch({
      py <- reticulate::import("pyfaidx")
      # print(py)
    }, error = function(e) {
      cat("Error: pyfaidx is not available.\n")
    })

    # sort gtf first
    gtf <- rtracklayer::import.gff(gtf_file,format = "gtf")

    # sort
    sorted_gtf <- gtf |> data.frame() |>
      dplyr::filter(type %in% type) |>
      dplyr::arrange(seqnames,gene_name,gene_id,transcript_id,start,end)

    # output
    output_name = paste(gtf_file,".sorted.gtf",sep = "")
    rtracklayer::export.gff(sorted_gtf,
                            con = output_name,
                            format = "gtf")

    # run code
    pyscript.path = system.file("extdata", "getSeq.py", package = "riboTransVis")
    reticulate::source_python(pyscript.path)
    reticulate::py$py_extractSequence(gtf_file = output_name,
                                      genome_file = genome_file,
                                      transcript_id = reticulate::r_to_py(transcript_id),
                                      new_id = reticulate::r_to_py(new_id),
                                      type = type,
                                      out_file = out_file)
  }
}






#' Prepare Transcript FASTA File for Analysis
#'
#' This function processes a transcript FASTA file by renaming the sequence headers
#' with a combination of transcript IDs and gene symbols for downstream analysis.
#' The modified FASTA file is generated as output and can be used in workflows such as
#' ribosome profiling or RNA-seq analysis.
#'
#' @param transcript_fa Path to the input FASTA file containing transcript sequences.
#'   The sequence headers must contain transcript IDs and gene symbols.
#' @param output_path Path to the directory where the modified FASTA file will be saved.
#'   The default is the current working directory (`"./"`).
#'
#' @return The function writes a modified FASTA file called `transcript.fa` in the
#'   specified `output_path`. This file contains sequences with updated headers
#'   combining transcript IDs and gene symbols. Returns nothing to the R environment.
#'
#' @details
#' This function reads a transcript FASTA file, extracts transcript IDs and gene symbols
#' from the sequence headers, and renames the headers in the format `transcriptID|geneSymbol`.
#' The resulting FASTA file (`transcript.fa`) is saved for downstream analysis.
#'
#' @importFrom Biostrings readDNAStringSet writeXStringSet
#' @export
prepare_transcipt_file <- function(transcript_fa = NULL,output_path = "./"){
  fa <- Biostrings::readDNAStringSet(transcript_fa)

  # assign new symbol
  tidtmp <- sapply(strsplit(names(fa),split = " "),"[",1)
  tid <- sapply(strsplit(tidtmp,split = "\\."),"[",1)
  tmp <- sapply(strsplit(names(fa),split = "gene_symbol:"),"[",2)
  gname <- sapply(strsplit(tmp,split = " "),"[",1)

  newid <- paste(tid,gname,sep = "|")

  names(fa) <- newid

  # output
  Biostrings::writeXStringSet(fa,filepath = paste(output_path, "transcript.fa",sep = ""),
                              format = "fasta")

}
