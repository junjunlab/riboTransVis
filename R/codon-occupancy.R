#' Codon Occupancy Plot
#'
#' This function calculates the codon occupancy for ribosome profiling data and generates a bar plot.
#'
#' @param object An object of class \code{ribotrans}, containing ribosome profiling data.
#' @param merge_rep Logical. Whether to merge replicate samples by \code{sample_group}. Default is \code{FALSE}.
#' @param cds_fa Character. The path to a FASTA file containing CDS sequences.
#' @param do_offset_correct Logical. If `TRUE`, performs **offset correction**
#' using `do_offset_correction()`. **Default**: `FALSE`.
#' @param position_shift Integer defining how much to adjust **ribosome footprint positions**
#' during offset correction. **Default**: `0`.
#' @param exclude_length A numeric vector of length two specifying the number of nucleotides
#' to exclude from the start and end of the CDS to avoid boundary effects. Default is `c(90, 90)`.
#' @param min_counts A numeric value indicating the minimum read counts required for a transcript
#' to be included in the analysis. Default is `32`.
#' @param plot_abbreviation Logical. If \code{TRUE}, uses amino acid abbreviations instead of codons as x-axis labels. Default is \code{FALSE}.
#' @param facet Logical. If \code{TRUE}, the plot is faceted by amino acid groups. Default is \code{TRUE}.
#' @param return_data Logical. If \code{TRUE}, returns a data frame with codon occupancy data instead of generating a plot. Default is \code{FALSE}.
#' @param ... Additional arguments (currently unused).
#'
#' @return
#' If \code{return_data = FALSE}, returns a ggplot2 object (bar plot of codon occupancy).
#' If \code{return_data = TRUE}, returns a data frame with the following columns:
#' \itemize{
#'   \item \code{sample} - Sample identifier.
#'   \item \code{codon_seq} - Codon sequence.
#'   \item \code{AminoAcid} - Full name of the amino acid.
#'   \item \code{Abbreviation3} - Three-letter amino acid abbreviation.
#'   \item \code{Abbreviation1} - One-letter amino acid abbreviation.
#'   \item \code{reloccup} - Relative codon occupancy.
#' }
#'
#' @details
#' Codon occupancy is calculated based on normalized ribosome profiling read coverage. The function:
#' \enumerate{
#'   \item Filters out invalid CDS data.
#'   \item Computes average read counts per codon.
#'   \item Loads CDS sequences and extracts codon-level information.
#'   \item Merges ribosome occupancy data with codon identity.
#'   \item Generates a relative codon occupancy profile.
#'   \item Optionally, visualizes codon occupancy as a bar plot.
#' }
#'
#' @examples
#' \dontrun{
#' codon_ocuupancy_plot(object = my_ribotrans,
#'                      cds_fa = "cds.fa",
#'                      plot_abbreviation = TRUE,
#'                      facet = TRUE)
#' }
#'
#' @importFrom Biostrings readDNAStringSet
#' @importFrom fastplyr f_filter f_group_by f_summarise f_left_join
#' @import ggplot2
#' @import dplyr
#' @export
setGeneric("codon_occupancy_plot",function(object,...) standardGeneric("codon_occupancy_plot"))


#' @rdname codon_occupancy_plot
#' @export
setMethod("codon_occupancy_plot",
          signature(object = "ribotrans"),
          function(object,
                   merge_rep = FALSE,
                   cds_fa = NULL,
                   do_offset_correct = FALSE,
                   position_shift = 0,
                   exclude_length = c(90,90),
                   min_counts = 32,
                   plot_abbreviation = FALSE,
                   facet = TRUE,
                   return_data = FALSE){

            # whether do reads offset correction
            if(do_offset_correct == TRUE){
              sry <- do_offset_correction(object = object,shift = position_shift)
            }else{
              sry <- object@summary_info
            }

            sry <- sry %>% fastplyr::f_filter(mstart != 0 | mstop != 0)


            # expect reads per position
            avg.ct <- sry %>%
              dplyr::mutate(relst = pos - mstart, relsp = pos - mstop) %>%
              fastplyr::f_filter(relst > exclude_length[1] & relsp < -exclude_length[2]) %>%
              dplyr::mutate(cdslen = mstop - mstart + 1) %>%
              dplyr::group_by(sample,rname,cdslen) %>%
              dplyr::summarise(counts = sum(count)) %>%
              dplyr::mutate(avg_ct = counts/cdslen) %>%
              dplyr::select(sample,rname,counts,avg_ct)

            # average reads
            pltdf <- sry %>%
              dplyr::inner_join(y = avg.ct,by = c("sample", "rname")) %>%
              # filter low counts
              fastplyr::f_filter(counts > min_counts) %>%
              dplyr::mutate(rel = pos - mstart, norm = count/avg_ct) %>%
              fastplyr::f_filter(rel >= 0 & rel <= (mstop - mstart + 1)) %>%
              fastplyr::f_group_by(sample,sample_group,rname,rel) %>%
              fastplyr::f_summarise(normsm = sum(norm)) %>%
              # codon position
              dplyr::mutate(rel = (rel %/% 3) + 1) %>%
              fastplyr::f_group_by(sample,sample_group,rname,rel) %>%
              fastplyr::f_summarise(value = mean(normsm))

            # ==================================================================
            # filter cds sequence

            # load cds fasta
            cds <- Biostrings::readDNAStringSet(cds_fa)

            # remove transcript contains Ns in sequence
            valid_cds <- cds[!grepl("[^ACGTacgt]", as.character(cds))]

            ids_retain <- intersect(unique(sry$rname),names(valid_cds))
            valid_cds <- valid_cds[ids_retain]

            # Filter sequences whose length is a multiple of 3
            valid_cds <- valid_cds[Biostrings::width(valid_cds) %% 3 == 0]

            # loop to extract codon seqeunce
            lapply(seq_along(valid_cds),function(x){
              tmp <- valid_cds[x]
              cdslen <- Biostrings::width(tmp)
              interval <- seq(1, cdslen, by = 3)

              seqs <- stringr::str_sub_all(tmp,start = interval, end = interval + 2)[[1]]

              data.frame(rname = names(tmp),rel = 1:(cdslen/3),codon_seq = seqs)
            }) %>% do.call("rbind",.) %>% data.frame() -> codon_info


            # ==================================================================


            # merge with occupancy
            pltdf2 <- pltdf %>%
              fastplyr::f_inner_join(y = codon_info,by = c("rname", "rel")) %>%
              fastplyr::f_group_by(sample,sample_group,codon_seq) %>%
              fastplyr::f_summarise(occup = sum(value),freq = dplyr::n()) %>%
              dplyr::mutate(reloccup = occup/freq) %>% na.omit()


            # whether aggregate replicates
            if(merge_rep == TRUE){
              pltdf2 <- pltdf2 %>%
                fastplyr::f_group_by(sample_group,codon_seq) %>%
                fastplyr::f_summarise(reloccup = mean(reloccup)) %>%
                dplyr::rename(sample = sample_group)
            }

            # amino acid annotation
            aa_info <- get_aa_table()

            # annotatate with amino acid info
            pltdf3 <- pltdf2 %>%
              fastplyr::f_inner_join(y = aa_info,by = c("codon_seq" = "Codon")) %>%
              dplyr::mutate(codon = paste(codon_seq, Abbreviation1,sep = " | "),
                            abbrev = paste(Abbreviation3, Abbreviation1,sep = " | "),)

            # check plot type
            if(plot_abbreviation == TRUE){
              pdf <- pltdf3 %>%
                dplyr::group_by(sample,abbrev,AminoAcid) %>%
                fastplyr::f_summarise(reloccup = sum(reloccup))

              lyr <- geom_col(aes(x = abbrev,y = reloccup, fill = sample),
                              position = position_dodge2())
              angle <- 90
            }else{
              pdf <- pltdf3
              lyr <- geom_col(aes(x = codon,y = reloccup, fill = sample),
                              position = position_dodge2())
              angle <- 90
            }

            # check facet
            if(facet == TRUE){
              facetlyr <- facet_wrap(~AminoAcid,scales = "free")
              angle <- 0
            }else{
              facetlyr <- NULL
            }
            # ==================================================================
            # plot
            # ==================================================================
            p <-
              ggplot(pdf) +
              lyr +
              theme_bw() +
              theme(axis.text.x = element_text(angle = angle,vjust = 0.5),
                    panel.grid = element_blank(),
                    strip.text = element_text(face = "bold"),
                    axis.text = element_text(colour = "black")) +
              xlab("Codons (Amino acids)") +
              ylab("Codon occupancy") +
              facetlyr

            # return
            if(return_data == FALSE){
              return(p)
            }else{
              return(pdf)
            }
          }
)






# ==============================================================================
# codon occupancy scatter plot
# ==============================================================================

#' Codon or Amino Acid Occupancy Scatter Plot
#'
#' Generate a scatter plot comparing codon or amino acid occupancy between two samples.
#' Supports optional labeling of specific codons/amino acids, with identity line and correlation statistics.
#'
#' @param codon_data A \code{data.frame} or \code{tibble} containing occupancy data. Should include columns \code{sample}, \code{codon_seq} (for codon analysis), \code{Abbreviation1} (for amino acid analysis), and \code{reloccup}.
#' @param x Character. Sample name to be shown on the x-axis (must match one of the sample names in \code{codon_data}).
#' @param y Character. Sample name to be shown on the y-axis (must match one of the sample names in \code{codon_data}).
#' @param type Character. Type of sequence to analyze: either \code{"codon"} or \code{"amino"}. Default is \code{"codon"}.
#' @param codon_labels Character vector. Optional. A vector of codon triplets or one-letter amino acid codes to label in the plot.
#' @param color Character. Color code for points (default: \code{"#FF6600"}).
#' @param size Numeric. Size of points in the scatter plot (default: 2).
#'
#' @return A \code{ggplot2} object of the scatter plot, showing codon or amino acid occupancy in two conditions.
#'
#' @details
#' Function workflow:
#' \itemize{
#'   \item For type = "codon": reshapes \code{codon_data} to wide format and plots codon-level scatter plot.
#'   \item For type = "amino": averages occupancy per amino acid and plots amino acid-level scatter plot.
#'   \item Optionally labels selected codons or amino acids.
#'   \item Includes reference y = x dashed line and Pearson correlation coefficient (via \code{ggpubr::stat_cor}).
#' }
#'
#'
#' @examples
#' \dontrun{
#' # For codon-level comparison
#' codon_scatter_plot(codon_data = df, x = "wt", y = "mut", type = "codon",
#'                    codon_labels = c("AAA", "GGG"))
#'
#' # For amino acid comparison
#' codon_scatter_plot(codon_data = df, x = "ctrl", y = "treat", type = "amino")
#' }
#'
#' @export
codon_scatter_plot <- function(codon_data = NULL,
                               x = NULL, y = NULL,
                               type = c("codon","amino"),
                               codon_labels = NULL,
                               color = "#FF6600", size = 2){
  type <- match.arg(type,choices = c("codon","amino"))
  # ============================================================================
  if(type == "codon"){
    ccdf_wide <- codon_data %>%
      dplyr::select(sample, codon_seq, reloccup) %>%
      tidyr::pivot_wider(names_from = sample, values_from = reloccup)

    ccdf_wide[is.na(ccdf_wide)] <- 0

    # mark label
    if(!is.null(codon_labels)){
      mk <- ggrepel::geom_text_repel(data = subset(ccdf_wide, codon_seq %in% codon_labels),
                                     aes(label = codon_seq),max.overlaps = 200)
    }else{
      mk <- ggrepel::geom_text_repel(aes(label = codon_seq),max.overlaps = 200)
    }

    # scatter plot
    p <-
      ggplot(ccdf_wide,
             aes(x = .data[[x]],y = .data[[y]])) +
      geom_abline(slope = 1,intercept = 0,lty = "dashed") +
      geom_point(size = size,color = color) +
      mk +
      ggtitle("Codon occupancy") +
      theme_bw(base_size = 12) +
      theme(panel.grid = element_blank(),
            plot.title = element_text(hjust = 0.5),
            axis.text = element_text(colour = "black")) +
      xlim(0,round(max(ccdf$reloccup),digits = 0)) + ylim(0,round(max(ccdf$reloccup),digits = 0)) +
      coord_equal() +
      ggpubr::stat_cor()
  }else{
    # ==========================================================================
    # aa
    ccdf_wide <- codon_data %>%
      dplyr::group_by(sample,Abbreviation1) %>%
      dplyr::summarise(reloccup = mean(reloccup)) %>%
      tidyr::pivot_wider(names_from = sample,values_from = reloccup)

    ccdf_wide[is.na(ccdf_wide)] <- 0

    if(!is.null(codon_labels)){
      mk <- ggrepel::geom_text_repel(data = subset(ccdf_wide, Abbreviation1 %in% codon_labels),
                                     aes(label = Abbreviation1),max.overlaps = 200)
    }else{
      mk <- ggrepel::geom_text_repel(aes(label = Abbreviation1),max.overlaps = 200)
    }

    # 5x5
    p <-
      ggplot(ccdf_wide,
             aes(x = .data[[x]],y = .data[[y]])) +
      geom_abline(slope = 1,intercept = 0,lty = "dashed") +
      geom_point(size = size,color = color) +
      mk +
      ggtitle("Amino acid occupancy") +
      theme_bw(base_size = 12) +
      theme(panel.grid = element_blank(),
            plot.title = element_text(hjust = 0.5),
            axis.text = element_text(colour = "black")) +
      xlim(0,round(max(ccdf_wide[,-1]),digits = 0)) + ylim(0,round(max(ccdf_wide[,-1]),digits = 0)) +
      coord_equal() +
      ggpubr::stat_cor()
  }

  return(p)
}
