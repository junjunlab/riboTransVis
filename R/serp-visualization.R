# ==============================================================================
# enrichment plot for single gene
# ==============================================================================

#' Plot Enrichment Ratios from Selective Ribosome Profiling
#'
#' @description
#' Generic function for plotting enrichment ratios from selective ribosome profiling experiments.
#'
#' @param object An S4 object containing enrichment data.
#' @param ... Additional arguments passed to specific method implementations.
#'
#' @return A plot visualizing the enrichment data.
#'
#' @details
#' This is a generic function that creates visualizations of enrichment ratios between
#' immunoprecipitated (IP) and total ribosome profiling samples. The specific implementation
#' depends on the class of the input object.
#'
#' @examples
#' \dontrun{
#' # Process a serp object
#' serp_obj <- get_occupancy(serp_obj, serp_exp = "total")
#' serp_obj <- get_occupancy(serp_obj, serp_exp = "ip")
#' serp_obj <- get_enrichment(serp_obj)
#'
#' # Create enrichment plot
#' enrichment_plot(serp_obj)
#' }
#'
#' @export
setGeneric("enrichment_plot",function(object,...) standardGeneric("enrichment_plot"))






#' @describeIn enrichment_plot Plot enrichment ratios for SERP objects
#'
#' @param cds_line Vector. Color and line width for the CDS indicator line. Default: c("#CC6600", 4).
#' @param line_col Character. Color of the enrichment line. Default: "#377CC0" (blue).
#' @param shadow_alpha Numeric vector. Alpha range for confidence interval shading. Default: c(0.1, 0.3).
#' @param enrich_threshold Numeric. Threshold line for enrichment ratio. Default: 1.
#' @param new_signal_range Logical. Whether to display signal range annotation. Default: FALSE.
#' @param range_x Numeric. X position of signal range annotation (normalized). Default: 0.99.
#' @param range_y Numeric. Y position of signal range annotation (normalized). Default: 0.98.
#' @param range_size Numeric. Font size of signal range annotation. Default: 4.
#' @param range_digit Numeric. Number of decimal places for signal range annotation. Default: 1.
#' @param nrow Numeric. Number of rows when arranging multiple plots. Default: NULL.
#' @param ncol Numeric. Number of columns when arranging multiple plots. Default: NULL.
#'
#' @details
#' For objects of class 'serp', this method creates a visualization of the enrichment ratio
#' between IP and total ribosome profiling samples. The method requires enrichment data to be
#' present in the 'enriched_ratio' slot of the object, which can be calculated using the
#' `get_enrichment` method.
#'
#' The plot includes:
#' 1. A line representing the mean enrichment ratio at each position
#' 2. Shaded areas representing confidence intervals, with opacity reflecting precision
#' 3. A dashed horizontal line at the specified enrichment threshold (typically 1)
#' 4. A bar at the bottom indicating the coding sequence (CDS) region
#' 5. Optional signal range annotation
#'
#' The x-axis is scaled to represent codon/amino acid positions rather than nucleotide positions.
#'
#' @return A ggplot object or a combined plot grid if multiple transcripts are present.
#'
#' @examples
#' \dontrun{
#' # Process a serp object
#' serp_obj <- get_occupancy(serp_obj, serp_exp = "total")
#' serp_obj <- get_occupancy(serp_obj, serp_exp = "ip")
#' serp_obj <- get_enrichment(serp_obj)
#'
#' # Basic enrichment plot
#' enrichment_plot(serp_obj)
#'
#' # Customized enrichment plot
#' enrichment_plot(serp_obj,
#'                line_col = "darkred",
#'                shadow_alpha = c(0.05, 0.2),
#'                enrich_threshold = 1.5,
#'                new_signal_range = TRUE)
#'
#' # Multiple transcripts in a grid
#' enrichment_plot(serp_obj, nrow = 2, ncol = 3)
#' }
#'
#' @importFrom ggplot2 ggplot geom_rect geom_line geom_hline scale_alpha_continuous facet_grid theme element_blank element_text xlab ylab
#' @importFrom ggside geom_xsidesegment scale_xsidey_continuous ggside
setMethod("enrichment_plot",
          signature(object = "serp"),
          function(object,
                   cds_line = c("#CC6600",4),
                   line_col = "#377CC0",
                   shadow_alpha = c(0.1,0.3),
                   enrich_threshold = 1,
                   new_signal_range = F,
                   range_x = 0.99,
                   range_y = 0.98,
                   range_size = 4,
                   range_digit = 1,
                   nrow = NULL, ncol = NULL){
            # check data
            if(nrow(object@enriched_ratio) == 0){
              stop("Please run `get_enrichment` first!")
            }

            plotdf <- object@enriched_ratio

            # check data
            feature <- subset(object@features, gene == object@gene_name)

            # window size
            bin <- plotdf$window_size[1]

            tid <- feature$idnew

            # loop for each transcript_id
            # x = 1
            lapply(seq_along(tid),function(x){
              tmp_feature <- subset(feature, idnew == tid[x])
              tmp <- subset(plotdf, rname == tid[x])

              # ================================================================
              # whether add new signal range
              if(new_signal_range == TRUE){
                axis.text.y = element_blank()
                axis.ticks.y = element_blank()

                range <- paste("[0-",round(max(tmp$upper),digits = range_digit),"]",sep = "")

                if (requireNamespace("ggpp", quietly = TRUE)) {
                  range_label <- ggpp::annotate(geom = "text_npc",
                                                npcx = range_x,npcy = range_y,
                                                label = range,
                                                size = range_size)
                } else {
                  warning("Package 'ggpp' is needed for this function to work.")
                }

              }else{
                axis.text.y = NULL
                axis.ticks.y = NULL
                range_label <- NULL
              }

              # ================================================================
              # plot
              ggplot(tmp) +
                geom_rect(aes(xmin = (pos+bin/2)/3-0.5,
                              xmax = (pos+bin/2)/3+0.5,
                              ymin = lower,
                              ymax = upper,
                              alpha = 1/((1/win_total+1/win_ip)*bin) ),
                          fill = line_col,show.legend = F) +
                geom_line(aes(x = (pos+bin/2)/3, y = enrich), col = line_col) +
                geom_hline(yintercept = enrich_threshold,lty = "dashed",color = "black") +
                range_label +
                # scale_y_continuous(trans = "log2", oob = scales::squish, expand = c(0,0)) +
                scale_alpha_continuous(trans = "sqrt", name = "prec", range = shadow_alpha) +
                facet_grid(sample~rname,switch = "y") +
                theme(panel.grid = element_blank(),
                      axis.text = element_text(colour = "black"),
                      strip.text.y.left = element_text(angle = 0),
                      strip.background = element_blank(),
                      strip.text = element_text(face = "bold"),
                      strip.placement = "outside",
                      axis.text.y = axis.text.y,
                      axis.ticks.y = axis.ticks.y,
                      ggside.panel.background = element_blank(),
                      ggside.panel.border = element_blank()) +
                xlab("Ribosome position\n(codons / amino acids)") +
                ylab("Mean enrichment\n(IP / total)") +
                ggside::geom_xsidesegment(data = data.frame(x = 1,xend = tmp_feature$cds/3),
                                          aes(x = x,xend = xend,
                                              y = 0.5,yend = 0.5),
                                          inherit.aes = F,
                                          linewidth = as.numeric(cds_line[2]),color = cds_line[1]) +
                ggside::scale_xsidey_continuous(breaks = NULL) +
                # scale_x_continuous(expand = c(0,0)) +
                ggside::ggside(collapse = "x")


            }) -> plist

            # combine plot list
            if (requireNamespace("cowplot", quietly = TRUE)) {
              cowplot::plot_grid(plotlist = plist,nrow = nrow,ncol = ncol)
            } else {
              warning("Package 'cowplot' is needed for this function to work.")
            }
          }
)
