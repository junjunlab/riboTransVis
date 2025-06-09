# ==============================================================================
# define digestion site motif plot
# ==============================================================================

#' Visualize Nucleotide Motifs Around Digestion Sites
#'
#' This method generates a sequence motif plot of digestion sites for ribosome profiling data.
#' The motif is derived from sequences flanking the 5' or 3' end of aligned reads in transcriptomic regions.
#' The function supports frame-specific motif analysis, customizable motif size, and layout adjustment.
#'
#' @param object A \code{ribotrans} object containing Ribo-seq summary information.
#' @param transcript_fa A path to a FASTA file of reference transcript sequences. This file is used to extract nucleotide sequences around the digestion site positions.
#' @param read_length A numeric vector of length 2 specifying the minimum and maximum read length to consider (default: \code{c(20, 35)}).
#' @param digest_site_range A numeric vector of length 2 specifying the number of nucleotides to extract upstream (5') and downstream (3') of the digestion site (default: \code{c(8, 8)}).
#' @param type A character string specifying which digestion end to use for motif extraction; either \code{"end5"} or \code{"end3"}.
#' @param show_frame Logical. Whether to split motif plots by ribosomal frame (0, 1, 2). Default is \code{FALSE}.
#' @param method A character string specifying the logo scaling method. Either \code{"prob"} for position probability matrices or \code{"bits"} for information content (bits). Passed to \code{ggseqlogo}.
#' @param merge_rep Logical. Whether to merge replicates based on \code{sample_group}. Default is \code{FALSE}.
#' @param return_data Logical. If \code{TRUE}, returns raw position frequency matrices (PFMs) instead of ggplot objects. Default is \code{FALSE}.
#' @param ... Additional arguments (currently unused).
#'
#'
#' @return Returns either:
#' \itemize{
#'   \item A nested list of ggplot objects (one per sample, optionally separated by frame), for logo visualization.
#'   \item A list of PFMs (if \code{return_data = TRUE}).
#' }
#'
#' @details
#' This function filters input Ribo-seq reads based on read length and alignment annotation (using \code{mstart} and \code{mstop}),
#' calculates digestion site coordinates depending on the selected \code{type}, extracts sequences from the corresponding transcriptome,
#' and computes a position frequency matrix (PFM) of nucleotides using \code{Biostrings::consensusMatrix}.
#' The logos are plotted using \code{ggseqlogo}, optionally highlighting the digestion site and breaking down by reading frame.
#'
#' @importFrom Rsamtools FaFile
#' @importFrom Biostrings getSeq consensusMatrix
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @importFrom dplyr filter mutate distinct rename
#' @importFrom fastplyr f_filter f_select
#'
#'
#' @examples
#' \dontrun{
#' # Load example ribotrans object and transcript FASTA file
#' digestion_site_plot(object = rb, transcript_fa = "trans.fa", type = "end5")
#'
#' # Show frame-specific motifs in bit score
#' digestion_site_plot(object = rb, transcript_fa = "trans.fa", show_frame = TRUE, method = "bits")
#' }
#'
#' @rdname digestion_site_plot
#' @export
setGeneric("digestion_site_plot",function(object,...) standardGeneric("digestion_site_plot"))





#' @rdname digestion_site_plot
#' @export
setMethod("digestion_site_plot",
          signature(object = "ribotrans"),
          function(object,
                   transcript_fa = NULL,
                   read_length = c(20,35),
                   digest_site_range = c(8,8),
                   type = c("end5", "end3"),
                   show_frame = FALSE,
                   method = c("prob","bits"),
                   merge_rep = FALSE,
                   return_data = FALSE){
            type <- match.arg(type,choices = c("end5", "end3"))
            method <- match.arg(method,choices = c("prob","bits"))
            # ==================================================================
            # filter data
            sry <- object@summary_info %>%
              fastplyr::f_filter(mstart != 0 | mstop != 0) %>%
              dplyr::filter(qwidth >= read_length[1] & qwidth <= read_length[2]) %>%
              dplyr::mutate(frame = (pos - mstart)%%3)

            if(merge_rep == TRUE){
              sry <- sry %>%
                fastplyr::f_select(sample_group , rname, pos, qwidth, translen, frame) %>%
                dplyr::distinct() %>%
                dplyr::rename(sample = sample_group)
            }else{
              ry <- sry %>%
                fastplyr::f_select(sample , rname, pos, qwidth, translen, frame) %>%
                dplyr::distinct()
            }


            # extend cut site
            if(object@assignment_mode == "end5"){
              sry <- sry %>%
                dplyr::mutate(pos5st = pos - digest_site_range[1], pos5sp = pos + digest_site_range[2],
                              pos3st = pos + qwidth - 1 - digest_site_range[1], pos3sp = pos + qwidth - 1 + digest_site_range[2]) %>%
                fastplyr::f_select(sample, rname, pos5st, pos5sp, pos3st, pos3sp, translen,frame)
            }else{
              sry <- sry %>%
                dplyr::mutate(pos3st = pos - digest_site_range[1], pos3sp = pos + digest_site_range[2],
                              pos5st = pos - qwidth + 1 - digest_site_range[1], pos5sp = pos - qwidth + 1 + digest_site_range[2]) %>%
                fastplyr::f_select(sample, rname, pos5st, pos5sp, pos3st, pos3sp, translen,frame)
            }

            sry.ft <- sry %>%
              fastplyr::f_filter(pos5st > 0 & pos3st > 0 & pos5sp > 0 & pos3sp > 0) %>%
              fastplyr::f_filter(pos5sp <= translen - 3 & pos3sp <= translen - 3)

            # ==================================================================
            # extract cut site sequence
            fasta_file <- Rsamtools::FaFile(transcript_fa)
            open(fasta_file)

            sp <- unique(sry.ft$sample)

            # x = 1
            lapply(seq_along(sp),function(x){
              tmp <- subset(sry.ft, sample == sp[x])

              # loop for frame
              fm <- c(0,1,2)

              # check whetehr show frame info
              if(show_frame == TRUE){
                lapply(seq_along(fm),function(f){
                  tmp2 <- subset(tmp, frame == fm[f])

                  if(type == "end5"){
                    gr <- GenomicRanges::GRanges(
                      seqnames = tmp2$rname,
                      ranges = IRanges::IRanges(start = tmp2$pos5st, end = tmp2$pos5sp)
                    )
                  }else{
                    gr <- GenomicRanges::GRanges(
                      seqnames = tmp2$rname,
                      ranges = IRanges::IRanges(start = tmp2$pos3st, end = tmp2$pos3sp)
                    )
                  }

                  cutseq <- Biostrings::getSeq(fasta_file, gr)

                  mat <- Biostrings::consensusMatrix(cutseq, as.prob = FALSE)
                  mat <- mat[intersect(c("A","G","C","T"), rownames(mat)), ]

                  return(mat)
                }) -> mat.list

                names(mat.list) <- c("frame0", "frame1", "frame2")
              }else{
                if(type == "end5"){
                  gr <- GenomicRanges::GRanges(
                    seqnames = tmp$rname,
                    ranges = IRanges::IRanges(start = tmp$pos5st, end = tmp$pos5sp)
                  )
                }else{
                  gr <- GenomicRanges::GRanges(
                    seqnames = tmp$rname,
                    ranges = IRanges::IRanges(start = tmp$pos3st, end = tmp$pos3sp)
                  )
                }

                cutseq <- Biostrings::getSeq(fasta_file, gr)

                mat <- Biostrings::consensusMatrix(cutseq, as.prob = FALSE)
                mat <- mat[intersect(c("A","G","C","T"), rownames(mat)), ]

                mat.list <- list(mat)
              }


              return(mat.list)
            }) -> pfm.list

            names(pfm.list) <- sp

            # ==================================================================
            # plot
            lapply(seq_along(pfm.list),function(x){
              tmp <- pfm.list[x]

              # j = 1
              lapply(seq_along(tmp[[1]]),function(j){
                pfm <- tmp[[1]][[j]]

                if (requireNamespace("ggseqlogo", quietly = TRUE)) {
                  if(show_frame == TRUE){
                    title <- paste(names(tmp),names(tmp[[1]][j]))
                  }else{
                    title <- names(tmp)
                  }

                  ggseqlogo::ggseqlogo(pfm, method = method) +
                    geom_rect(aes(xmin = digest_site_range[1] + 1 -0.5,xmax = digest_site_range[1] + 1 +0.5,ymin = -Inf,ymax = Inf),
                              fill = "orange",alpha = 0.25) +
                    scale_x_continuous(breaks = seq(1,digest_site_range[2] - -digest_site_range[1] + 1),
                                       labels = seq(-digest_site_range[1],digest_site_range[2])) +
                    ggtitle(title) +
                    theme(legend.position = "none",
                          axis.text = element_text(size = 12),
                          plot.title = element_text(hjust = 0.5),
                          axis.line.y = element_line(),
                          axis.line.x = element_line(),
                          axis.ticks.y = element_line(),
                          axis.ticks.x = element_line()) +
                    xlab(paste("Distance to ",type,"site"))
                } else {
                  warning("Package 'ggseqlogo' is needed for this function to work.")
                }

              }) -> plist

              return(plist)
            }) -> pfinal

            # return
            if(return_data == FALSE){
              return(pfinal)
            }else{
              return(pfm.list)
            }
          }

)
