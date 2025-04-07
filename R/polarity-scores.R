# ==============================================================================
# method
# ==============================================================================

#' Generate Polarity Score Density Plot
#'
#' @description Generates a polarity score density plot for `ribotrans` objects.
#' The plot visualizes the distribution of calculated polarity scores across different samples.
#'
#' @param object A `ribotrans` object containing sequencing and ribosome profiling data.
#' @param merge_rep Logical. Whether to merge replicate samples by \code{sample_group}. Default is \code{FALSE}.
#' @param do_offset_correct Logical. If `TRUE`, performs **offset correction**
#' using `do_offset_correction()`. **Default**: `FALSE`.
#' @param position_shift Integer defining how much to adjust **ribosome footprint positions**
#' during offset correction. **Default**: `0`.
#' @param exclude_length Numeric vector of length two, specifying the number of nucleotides
#' to exclude from the start and end of coding sequences. (Default: `c(15,15)`)
#' @param min_counts Numeric, minimum threshold for the total read count of a gene to be included in the analysis. (Default: `64`)
#' @param mark_line_col Character string defining the color of the vertical guiding lines marking key positions. (Default: `"grey50"`)
#' @param orf_width Numeric, linewidth for the Open Reading Frame (ORF) annotation. (Default: `5`)
#' @param orf_col Character string defining the color for the ORF annotation. (Default: `"grey"`)
#' @param return_data Logical, if `TRUE`, returns processed density data instead of a plot. (Default: `FALSE`)
#' @param ... Additional arguments (currently unused).
#'
#' @return
#' If `return_data = FALSE`, returns a `ggplot2` object visualizing the polarity score density by sample.
#' If `return_data = TRUE`, returns a `data.frame` containing calculated density values.
#'
#' @details
#' - Reads are normalized to **Reads Per Million (RPM)** mapped reads.
#' - Genes with read counts below `min_counts` are excluded from calculations.
#' - Density is estimated using **kernel density estimation (KDE)** with `density(..., adjust = 1)`.
#' - Vertical dashed lines represent key regions: `x = 0` (balanced polarity) and `max density peak` per sample.
#' - The ORF region is annotated as a reference.
#'
#' @examples
#' \dontrun{
#' # Generate polarity plot for a `ribotrans` object
#' polarity_plot(ribo_obj)
#'
#' # Retrieve the processed density data instead of a plot
#' data <- polarity_plot(ribo_obj, return_data = TRUE)
#' head(data)
#' }
#'
#' @import ggplot2
#' @import dplyr
#' @importFrom fastplyr f_filter f_select f_group_by f_summarise f_left_join
#' @import ggside
#' @importFrom purrr map_df
#' @importFrom stats density
#' @export
setGeneric("polarity_plot",function(object,...) standardGeneric("polarity_plot"))



#' @rdname polarity_plot
#' @export
setMethod("polarity_plot",
          signature(object = "ribotrans"),
          function(object,
                   merge_rep = FALSE,
                   do_offset_correct = FALSE,
                   position_shift = 0,
                   exclude_length = c(15,15),
                   min_counts = 64,
                   mark_line_col = "grey50",
                   orf_width = 5,
                   orf_col = "grey",
                   return_data = FALSE){
            features <- object@features

            # whether do reads offset correction
            if(do_offset_correct == TRUE){
              sry <- do_offset_correction(object = object,shift = position_shift)
            }else{
              sry <- object@summary_info
            }

            # exclude cds first and last nts
            sry <- sry %>%
              fastplyr::f_filter(mstart != 0 | mstop != 0) %>%
              fastplyr::f_select(-qwidth,-translen) %>%
              dplyr::mutate(relst = pos - mstart,
                            relsp = pos - mstop) %>%
              fastplyr::f_filter(relst > exclude_length[1] & relsp < -exclude_length[2])

            # rpm normalization
            lib <- subset(object@library, type == "ribo") %>%
              dplyr::select(sample, mappped_reads)

            sry <- sry %>%
              dplyr::left_join(y = lib,by = "sample") %>%
              dplyr::mutate(rpm = (count/mappped_reads)*10^6)

            # filter low counts
            density.tt <- sry %>%
              fastplyr::f_group_by(sample,rname) %>%
              fastplyr::f_summarise(count.tt = sum(count),
                                    rpm.tt = sum(rpm)) %>%
              fastplyr::f_filter(count.tt > min_counts)

            # calculate polarity scores
            pl_df <- density.tt |>
              fastplyr::f_left_join(y = sry,by = c("sample","rname")) |>
              dplyr::mutate(wi = (2*relst - (mstop - mstart + 1))/(mstop - mstart - 1),
                            pi = rpm*wi/rpm.tt) |>
              fastplyr::f_group_by(sample,sample_group,rname) |>
              fastplyr::f_summarise(sum_pi = sum(pi))

            # whether aggregate replicates
            if(merge_rep == TRUE){
              pl_df <- pl_df %>%
                fastplyr::f_group_by(sample_group,rname) %>%
                fastplyr::f_summarise(sum_pi = mean(sum_pi)) %>%
                dplyr::rename(sample = sample_group)
            }

            # calculate density

            sp <- unique(pl_df$sample)

            purrr::map_df(seq_along(sp),function(x){
              tmp <- subset(pl_df, sample == sp[x])

              density_data <- stats::density(tmp$sum_pi, na.rm = TRUE, adjust = 1)

              # bw <- density_data$bw
              df <- data.frame(x = density_data$x,
                               y = density_data$y * nrow(tmp))
              df$sample <- sp[x]

              return(df)
            }) -> density.df

            # filter max position
            ymax <- density.df %>%
              dplyr::slice_max(order_by = y,by = sample)

            # ==================================================================
            # plot
            # ==================================================================
            struc <- data.frame(x = -1,xend = 1,y = 0.5,yend = 0.5)
            struc.anno <- data.frame(x = 0,y = 0.5,label = c("ORF"))

            p <-
              ggplot(density.df) +
              geom_path(aes(x = x,y = y)) +
              geom_vline(xintercept = 0,lty = "dashed",color = mark_line_col) +
              geom_vline(data = ymax,aes(xintercept = x),lty = "dashed",color = mark_line_col) +
              facet_wrap(~sample) +
              xlab("Polarity score") +
              ylab("Number of genes") +
              scale_x_continuous(breaks = c(-1,-0.5,0,0.5,1),
                                 limits = c(-1,1),
                                 labels = c("-1","-0.5","0","0.5","1")) +
              # annotation
              geom_xsidesegment(data = struc,aes(x = x,xend = xend,y = y,yend = yend),
                                inherit.aes = F,
                                color = orf_col,
                                linewidth = orf_width) +
              geom_xsidetext(data = struc.anno,aes(x = x,y = y,label = label)) +
              theme(ggside.panel.background = element_blank(),
                    ggside.panel.border = element_blank(),
                    panel.grid = element_blank(),
                    axis.text = element_text(colour = "black"),
                    strip.text = element_text(colour = "black",size = rel(1))) +
              ggside::scale_xsidey_continuous(breaks = NULL)

            # return
            if(return_data == FALSE){
              return(p)
            }else{
              return(density.df)
            }
          }
)
