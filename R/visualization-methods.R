# ==============================================================================
# visualization
# ==============================================================================


#' Plot RNA-seq and Ribosome Occupancy along Transcripts
#'
#' @description
#' Visualizes ribosome footprint occupancy or RNA coverage across transcript coordinates for a given gene using data stored in a \code{ribotrans} object. This function supports various data layers (e.g., IP, RNA, scaled ribo), replicate visualization, structural annotations, and codon/nucleotide views.
#'
#' @param object A \code{ribotrans} object. Must have properly populated occupancy slots (e.g., \code{ribo_occupancy}, \code{RNA_coverage}, etc.) and transcript \code{features}.
#' @param selected_id Character vector. A list of transcript IDs to plot. If \code{NULL} (default), all transcripts for the current gene will be plotted.
#' @param type Character. The type of signal to plot: one of \code{"ribo"}, \code{"rna"}, \code{"ribo_rna"}, or \code{"scaled_ribo"}. Default: \code{"ribo"}.
#' @param layer Layer geometry: \code{"line"} (default) for \code{geom_path()} curves, or \code{"col"} for \code{geom_col()} bar charts.
#' @param col_alpha Numeric in c(0,1), transparency for the bars when
#' \code{layer = "col"}. Default: \code{1}.
#' @param sample_order Character vector. Optional. The desired order of samples in the facet display.
#' @param facet_layer A ggplot2 facet object (e.g., \code{facet_grid()}) for organizing sample/transcript panels. Default: \code{facet_grid(sample ~ rname, switch = "y")}.
#' @param sep_mer_sample Logical. If \code{TRUE}, merged sample display is combined with individual sample facets. Default: \code{FALSE}.
#' @param new_signal_range Logical. Whether to annotate each panel with the signal range. Requires package \pkg{ggpp}. Default: \code{TRUE}.
#' @param position_mode Character. Plot resolution: \code{"nt"} (nucleotide) or \code{"codon"} (amino acid/codon level). Default: \code{"nt"}.
#' @param range_x X position (normalized 0–1) for the signal range label. Default: 0.98.
#' @param range_y Y position (normalized 0–1) for the signal range label. Default: 0.98.
#' @param range_size Numeric. Font size for signal range label. Default: 4.
#' @param range_digit Integer. Digits to round signal range values for annotation. Default: 1.
#' @param scale_factor Numeric. Scaling factor for RNA signal when comparing ribo and RNA together (used only when \code{type = "ribo_rna"}). Default: 1.
#' @param utr_width Line width for 5'/3' UTR region annotations respectively. Default: \code{1} for UTR.
#' @param cds_width Line width for CDS region annotations respectively. Default: \code{3} for CDS.
#' @param utr_col Colors used to annotate UTR in \code{ggside} alignment strip. Default: \code{"grey"}.
#' @param cds_col Colors used to annotate CDS in \code{ggside} alignment strip. Default: \code{"grey30"}.
#' @param nrow Number of rows in the final layout when multiple transcripts are plotted. Passed to \code{cowplot::plot_grid()}.
#' @param ncol Number of columns in the final layout when multiple transcripts are plotted. Passed to \code{cowplot::plot_grid()}.
#'
#' @details
#' This method supports four types of data display:
#'
#' \describe{
#' \item{"ribo"}{Ribosome footprint data from \code{ribo_occupancy} slot}
#' \item{"rna"}{Total RNA coverage data from \code{RNA_coverage} slot}
#' \item{"ribo_rna"}{Overlaid ribo + scaled RNA signals}
#' \item{"scaled_ribo"}{Pre-standardized/normalized ribo signal from \code{scaled_occupancy}}
#' }
#'
#' Transcript structure data is pulled from the \code{features} slot using the gene of interest (\code{object@gene_name}). For each transcript, smoothed signal trajectories across positions are displayed, with optional CDS annotation and signal range label.
#'
#' If \code{sep_mer_sample = TRUE}, individual samples and merged signal (e.g., mean profile) will be plotted together. Structural annotations such as CDS length and UTRs are visualized using \pkg{ggside}'s \code{geom_xsidesegment()}.
#'
#' @return A \code{ggplot} object (consisting of multiple panels combined using \pkg{cowplot}); each panel represents one transcript. Plot can show either separate or merged replicate signals.
#'
#' @seealso \code{\link{get_occupancy}}, \code{\link{get_coverage}}, \code{\link{get_scaled_occupancy}}
#'
#' @import ggplot2
#' @importFrom ggside geom_xsidesegment ggside scale_xsidey_continuous
#' @importFrom dplyr filter mutate group_by summarise rename
#'
#'
#' @examples
#' \dontrun{
#' # Plot ribosome occupancy with default settings
#' trans_plot(ribotrans_obj)
#'
#' # Plot RNA signal in bar chart format at nucleotide resolution
#' trans_plot(ribotrans_obj, type = "rna", layer = "col", position_mode = "nt")
#'
#' # Overlay Ribo & RNA with scaling + CDS annotation at codon level
#' trans_plot(ribotrans_obj, type = "ribo_rna", position_mode = "codon", scale_factor = 0.3)
#' }
#'
#' @export
setMethod("trans_plot",
          signature(object = "ribotrans"),
          function(object,
                   selected_id = NULL,
                   type = c("ribo","rna","ribo_rna","scaled_ribo"),
                   layer = c("line", "col"),
                   col_alpha = 1,
                   sample_order = NULL,
                   facet_layer = ggplot2::facet_grid(sample~rname,switch = "y"),
                   sep_mer_sample = FALSE,
                   new_signal_range = TRUE,
                   position_mode = c("nt", "codon"),
                   range_x = 0.98,
                   range_y = 0.98,
                   range_size = 4,
                   range_digit = 1,
                   scale_factor = 1,
                   utr_width = 1,
                   cds_width = 3,
                   utr_col = "grey",
                   cds_col = "grey30",
                   nrow = NULL,ncol = NULL){
            # check plot type
            type <- match.arg(type,choices = c("ribo","rna","ribo_rna","scaled_ribo"))
            layer <- match.arg(layer, choices = c("line", "col"))
            position_mode <- match.arg(position_mode, choices = c("nt", "codon"))

            if(type == "ribo"){
              # check data
              if(nrow(object@ribo_occupancy) == 0){
                stop("Please run `get_occupancy` first!")
              }

              pldf <- object@ribo_occupancy
              ylab <- "Ribsome footprint occupancy (RPM)"

              if(layer == "line"){
                player <- geom_path(aes(x = pos,y = smooth,color = sample))
              }else{
                player <- geom_col(aes(x = pos,y = smooth,fill = sample,color = sample),width = 1)
              }

              col <- NULL
              col1 <- NULL
            }else if(type == "rna"){
              # check data
              if(nrow(object@RNA_coverage) == 0){
                stop("Please run `get_coverage` first!")
              }

              pldf <- object@RNA_coverage
              ylab <- "RNA reads coverage (RPM)"

              if(layer == "line"){
                player <- geom_path(aes(x = pos,y = smooth,color = sample))
              }else{
                player <- geom_col(aes(x = pos,y = smooth,fill = sample,color = sample),width = 1)
              }

              col <- NULL
              col1 <- NULL
            }else if(type == "ribo_rna"){
              ylab <- "Ribsome footprint occupancy \n (RNA reads coverage)"

              # check data
              if(nrow(object@ribo_occupancy) == 0){
                stop("Please run `get_occupancy` first!")
              }

              # check data
              if(nrow(object@RNA_coverage) == 0){
                stop("Please run `get_coverage` first!")
              }

              ribo <- object@ribo_occupancy
              ribo$exp <- "ribo"
              rna <- object@RNA_coverage
              rna$exp <- "rna"
              rna$smooth <- rna$smooth*scale_factor

              pldf <- rbind(ribo, rna)
              pldf$exp <- factor(pldf$exp,levels = c("rna","ribo"))

              if(layer == "line"){
                player <- geom_path(aes(x = pos,y = smooth,color = exp))
                col <- scale_color_manual(values = c("ribo" = "red", "rna" = "grey"), name = "")
                col1 <- NULL
              }else{
                player <- geom_col(aes(x = pos,y = smooth,fill = exp,color = exp),width = 1, alpha = col_alpha)
                col <- scale_fill_manual(values = c("ribo" = "red", "rna" = "grey"), name = "")
                col1 <- scale_color_manual(values = c("ribo" = "red", "rna" = "grey"), name = "")
              }

            }else if(type == "scaled_ribo"){
              ylab <- "Scaled ribosome footprint occupancy"

              # check data
              if(nrow(object@scaled_occupancy) == 0){
                stop("Please run `get_scaled_occupancy` first!")
              }


              pldf <- object@scaled_occupancy

              if(layer == "line"){
                player <- geom_path(aes(x = pos,y = smooth,color = sample))
              }else{
                player <- geom_col(aes(x = pos,y = smooth,fill = sample,color = sample),width = 1)
              }

              col <- NULL
              col1 <- NULL
            }


            # filter genes
            tanno <- object@features %>%
              dplyr::filter(gene == object@gene_name) %>%
              dplyr::mutate(rname = idnew)

            if(is.null(selected_id)){
              tids <- tanno$idnew
            }else{
              tids <- selected_id
            }

            # ==================================================================
            # loop plot for each transcript
            # x = 1
            lapply(seq_along(tids),function(x){
              tmp.df <- subset(pldf, rname == tids[x])
              tanno.tmp <- subset(tanno, rname == tids[x])

              # check position_mode
              if(position_mode == "codon"){
                tmp.df <- tmp.df %>%
                  dplyr::filter(pos >= tanno.tmp$mstart & pos <= tanno.tmp$mstop) %>%
                  dplyr::mutate(pos = pos - tanno.tmp$mstart + 1) %>%
                  dplyr::mutate(codon = (pos - 1) %/% 3 + 1) %>%
                  dplyr::group_by(sample,rname,codon) %>%
                  dplyr::summarise(smooth = mean(smooth)) %>%
                  dplyr::rename(pos = codon)

                xlab <- "Position along transcript (codon/AA)"

                struc_layer1 <-
                  geom_xsidesegment(data = tanno.tmp,
                                    aes(x = 1,xend = tanno.tmp$cds/3,y = 0.5,yend = 0.5),
                                    linewidth = cds_width,color = cds_col)
                struc_layer2 <- NULL
              }else{
                xlab <- "Position along transcript (nt)"

                struc_layer1 <-
                  geom_xsidesegment(data = tanno.tmp,
                                    aes(x = 1,xend = translen,y = 0.5,yend = 0.5),
                                    linewidth = utr_width,color = utr_col)
                struc_layer2 <-
                  geom_xsidesegment(data = tanno.tmp,
                                    aes(x = mstart,xend = mstop,y = 0.5,yend = 0.5),
                                    linewidth = cds_width,color = cds_col)
              }

              # whether add new signal range
              if(new_signal_range == TRUE){
                axis.text.y = element_blank()
                axis.ticks.y = element_blank()

                range <- paste("[0-",round(max(tmp.df$smooth),digits = range_digit),"]",sep = "")

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

              # whether show merged sample and seperate sample together
              if(sep_mer_sample == TRUE){
                tmp.df2 <- tmp.df
                tmp.df2$sp <- "merged sample"

                tmp.df$sp <- tmp.df$sample
                tmp.df <- rbind(tmp.df2, tmp.df)

                facet_layer <- ggplot2::facet_grid(sp~rname,switch = "y")
              }

              # sample orders
              if(!is.null(sample_order)){
                if(sep_mer_sample == TRUE){
                  tmp.df$sp <- factor(tmp.df$sp ,levels = c(sample_order,"merged sample"))
                }else{
                  tmp.df$sample <- factor(tmp.df$sample ,levels = sample_order)
                }

              }

              # order
              if(type == "ribo_rna"){
                tmp.df <- tmp.df[order(tmp.df$exp=="ribo"), ]
              }


              # plot
              p <-
                ggplot(tmp.df) +
                player +
                range_label +
                theme_bw() +
                # facet_wrap(~sample,scales = "fixed",ncol = 1,switch = "y") +
                facet_layer +
                theme(panel.grid = element_blank(),
                      axis.text = element_text(colour = "black"),
                      strip.text = element_text(colour = "black",size = rel(1)),
                      strip.text.x = element_text(face = "bold",colour = "black",size = rel(1)),
                      strip.background = element_blank(),
                      strip.placement = "outside",
                      strip.text.y.left = element_text(angle = 0,hjust = 1),
                      axis.text.y = axis.text.y,
                      axis.ticks.y = axis.ticks.y,
                      ggside.panel.background = element_blank(),
                      ggside.panel.border = element_blank()) +
                # gene structure
                struc_layer1 + struc_layer2 +
                ggside::scale_xsidey_continuous(breaks = NULL) +
                xlab(xlab) +
                ylab(ylab) +
                ggside(collapse = "x") +
                col + col1

              return(p)
            }) -> plist

            # combine plot list
            if (requireNamespace("cowplot", quietly = TRUE)) {
              cowplot::plot_grid(plotlist = plist,nrow = nrow,ncol = ncol)
            } else {
              warning("Package 'cowplot' is needed for this function to work.")
            }


          }
)




# ==============================================================================
# function to visualize coverage and occupancy on genome
# ==============================================================================

#' Plot Transcript-level Ribosome and RNA Coverage with Gene Structure
#'
#' This function visualizes ribosome footprint occupancy and/or RNA read
#' coverage along transcripts for a given gene, together with the underlying
#' transcript structure (exons, CDS, UTRs). It returns a \code{ggplot} object
#' that can be further customized.
#'
#'
#' @param object A \code{ribotrans} object. Must contain the results of
#' \code{get_occupancy()} for ribosome data and/or \code{get_coverage()}
#' for RNA data, as well as the GTF annotation for the target gene stored
#' in \code{@gtf_data}.
#' @param selected_id Character vector of transcript IDs to include in the plot.
#' Default: \code{NULL}, in which case all transcripts for the gene are shown.
#' @param type Character; one of \code{"ribo"}, \code{"rna"} or
#' \code{"ribo_rna"}. Specifies whether to plot only ribosome occupancy,
#' only RNA coverage, or both. Default: \code{c("ribo", "rna", "ribo_rna")}
#' (i.e. \code{"ribo"}).
#' @param layer Character; one of \code{"col"} or \code{"line"}. Use bar plots
#' (\code{"col"}) or line plots (\code{"line"}) for the signal.
#' Default: \code{c("col", "line")} (i.e. \code{"col"}).
#' @param col_alpha Numeric in c(0,1), transparency for the bars when
#' \code{layer = "col"}. Default: \code{1}.
#' @param sample_order Optional character vector specifying the desired order
#' of samples in the facets. If \code{NULL}, the original order in the data
#' is retained. Default: \code{NULL}.
#' @param facet_layer A \code{ggplot2} facet specification (e.g.
#' \code{ggplot2::facet_grid}). Default:
#' \code{ggplot2::facet_grid(sample ~ rname, switch = "y")}.
#' @param sep_mer_sample Logical; if \code{TRUE}, displays both the merged
#' signal and individual samples side by side. Default: \code{FALSE}.
#' @param new_signal_range Logical; if \code{TRUE}, annotates a text label
#' showing the range of the plotted signal at position (\code{range_x},
#' \code{range_y}). Default: \code{FALSE}.
#' @param collapse_structure Logical; if \code{TRUE}, collapses all transcript
#' structures into a single line. Otherwise each transcript is plotted on
#' its own horizontal line. Default: \code{FALSE}.
#' @param range_x Numeric in c(0,1), x‐coordinate (in NPC units) of the range
#' annotation. Default: \code{0.9}.
#' @param range_y Numeric in c(0,1), y‐coordinate (in NPC units) of the range
#' annotation. Default: \code{0.9}.
#' @param range_size Numeric, text size of the range annotation. Default:
#' \code{4}.
#' @param range_digit Integer, number of digits to round the maximum signal
#' value for the label. Default: \code{1}.
#' @param scale_factor Numeric factor to scale the RNA coverage when
#' \code{type = "ribo_rna"}. Default: \code{1}.
#' @param utr_width Numeric, relative width of UTR segments (currently unused).
#' Default: \code{1}.
#' @param cds_width Numeric, relative width of CDS segments (currently unused).
#' Default: \code{3}.
#' @param background_line A length-2 vector \code{c(color, size)} controlling
#' the color and line width of the transcript backbone.
#' Default: \code{c("black", 0.25)}.
#' @param exon_line A length-2 vector \code{c(color, size)} controlling exon
#' segment styling. Default: \code{c("#003399", 1)}.
#' @param cds_line A length-2 vector \code{c(color, size)} controlling CDS
#' segment styling. Default: \code{c("#003399", 3)}.
#' @param structure_panel_height Numeric, relative height of the side panel
#' showing gene structure. Default: \code{0.2}.
#' @param ... Additional arguments (currently unused).
#'
#'
#' @return A \code{ggplot} object displaying the requested signal(s) along
#' the genomic coordinates, faceted by sample (and optionally merged
#' sample) and transcript region, with gene structure drawn at the side panel.
#'
#' @seealso
#' \code{\link{get_occupancy}}, \code{\link{get_coverage}} for preprocessing;
#' \code{ggplot2::facet_grid}, \code{ggside} for panel layout.
#'
#'
#' @export
setGeneric("genome_trans_plot",function(object,...) standardGeneric("genome_trans_plot"))




#' @rdname genome_trans_plot
#' @export
setMethod("genome_trans_plot",
          signature(object = "ribotrans"),
          function(object,
                   selected_id = NULL,
                   type = c("ribo","rna","ribo_rna"),
                   layer = c("col", "line"),
                   col_alpha = 1,
                   sample_order = NULL,
                   facet_layer = ggplot2::facet_grid(sample~rname,switch = "y"),
                   sep_mer_sample = FALSE,
                   new_signal_range = FALSE,
                   collapse_structure = FALSE,
                   range_x = 0.9,
                   range_y = 0.9,
                   range_size = 4,
                   range_digit = 1,
                   scale_factor = 1,
                   utr_width = 1,
                   cds_width = 3,
                   background_line = c("black",0.25),
                   exon_line = c("#003399",1),
                   cds_line = c("#003399",3),
                   structure_panel_height = 0.2){
            # check plot type
            type <- match.arg(type,choices = c("ribo","rna","ribo_rna"))
            layer <- match.arg(layer, choices = c("col", "line"))

            if(type == "ribo"){
              # check data
              if(nrow(object@ribo_occupancy) == 0){
                stop("Please run `get_occupancy` first!")
              }

              pldf <- object@ribo_occupancy
              ylab <- "Ribsome footprint occupancy (RPM)"

              if(layer == "line"){
                player <- geom_path(aes(x = pos,y = smooth,color = sample))
              }else{
                player <- geom_col(aes(x = pos,y = smooth,fill = sample,color = sample),width = 1)
              }

              col <- NULL
              col1 <- NULL
            }else if(type == "rna"){
              # check data
              if(nrow(object@RNA_coverage) == 0){
                stop("Please run `get_coverage` first!")
              }

              pldf <- object@RNA_coverage
              ylab <- "RNA reads coverage (RPM)"

              if(layer == "line"){
                player <- geom_path(aes(x = pos,y = smooth,color = sample))
              }else{
                player <- geom_col(aes(x = pos,y = smooth,fill = sample,color = sample),width = 1)
              }

              col <- NULL
              col1 <- NULL
            }else if(type == "ribo_rna"){
              ylab <- "Ribsome footprint occupancy \n (RNA reads coverage)"

              # check data
              if(nrow(object@ribo_occupancy) == 0){
                stop("Please run `get_occupancy` first!")
              }

              # check data
              if(nrow(object@RNA_coverage) == 0){
                stop("Please run `get_coverage` first!")
              }

              ribo <- object@ribo_occupancy
              ribo$exp <- "ribo"
              rna <- object@RNA_coverage
              rna$exp <- "rna"
              rna$smooth <- rna$smooth*scale_factor

              pldf <- rbind(ribo, rna)
              pldf$exp <- factor(pldf$exp,levels = c("rna","ribo"))

              if(layer == "line"){
                player <- geom_path(aes(x = pos,y = smooth,color = exp))
                col <- scale_color_manual(values = c("ribo" = "red", "rna" = "grey50"), name = "")
                col1 <- NULL
              }else{
                player <- geom_col(aes(x = pos,y = smooth,fill = exp,color = exp),width = 1,alpha = col_alpha)
                col <- scale_fill_manual(values = c("ribo" = "red", "rna" = "grey50"), name = "")
                col1 <- scale_color_manual(values = c("ribo" = "red", "rna" = "grey50"), name = "")
              }

            }


            # ==================================================================
            # structure layers
            # ==================================================================
            # add region structure
            gene_anno <- data.frame(object@gtf_data) %>%
              dplyr::filter(gene_name == object@gene_name)

            # check selected_id
            if(!is.null(selected_id)){
              gene_anno <- gene_anno %>%
                dplyr::filter(transcript_id %in% selected_id)
            }

            # add label
            chr <- ifelse(startsWith(as.character(gene_anno$seqnames[1]),"chr"),
                          as.character(gene_anno$seqnames[1]),
                          paste("chr",gene_anno$seqnames[1],sep = ""))

            pldf$rname <- paste(object@gene_name, " ",
                                chr,":",min(gene_anno$start),"-",max(gene_anno$end),sep = "")

            # ==================================================================
            # transcript background line
            sruc_trans_rg <- gene_anno %>%
              dplyr::filter(gene_name == object@gene_name & type == "transcript") %>%
              dplyr::mutate(group_id = dplyr::dense_rank(transcript_id))

            # check expand structure
            if(collapse_structure == TRUE){
              bg.layer <- geom_xsidesegment(data = sruc_trans_rg,
                                            aes(x = start,xend = end,y = 1,yend = 1),
                                            linewidth = as.numeric(background_line[2]),
                                            color = background_line[1])

              sturc_limit <- c(0,2)
            }else{
              bg.layer <- geom_xsidesegment(data = sruc_trans_rg,
                                            aes(x = start,xend = end,y = group_id,yend = group_id),
                                            linewidth = as.numeric(background_line[2]),
                                            color = background_line[1])

              sturc_limit <- c(0,nrow(sruc_trans_rg) + 1)
            }



            # ==================================================================
            # exons
            sruc_rg <- gene_anno %>%
              dplyr::filter(gene_name == object@gene_name & type == "exon") %>%
              dplyr::mutate(group_id = dplyr::dense_rank(transcript_id))

            # add exon structure
            if(collapse_structure == TRUE){
              exon.layer <- geom_xsidesegment(data = sruc_rg,
                                              aes(x = start,xend = end,y = 1,yend = 1),
                                              linewidth = as.numeric(exon_line[2]),
                                              color = exon_line[1])
            }else{
              exon.layer <- geom_xsidesegment(data = sruc_rg,
                                              aes(x = start,xend = end,y = group_id,yend = group_id),
                                              linewidth = as.numeric(exon_line[2]),
                                              color = exon_line[1])
            }


            # ================================================================
            # cds
            sruc_cds_rg <- gene_anno %>%
              dplyr::filter(gene_name == object@gene_name & type == "CDS") %>%
              dplyr::mutate(group_id = dplyr::dense_rank(transcript_id))

            # add cds structure
            if(nrow(sruc_cds_rg) > 0){
              if(collapse_structure == TRUE){
                cds.layer <- geom_xsidesegment(data = sruc_cds_rg,
                                               aes(x = start,xend = end,y = 1,yend = 1),
                                               linewidth = as.numeric(cds_line[2]),
                                               color = cds_line[1])
              }else{
                cds.layer <- geom_xsidesegment(data = sruc_cds_rg,
                                               aes(x = start,xend = end,y = group_id,yend = group_id),
                                               linewidth = as.numeric(cds_line[2]),
                                               color = cds_line[1])
              }

            }else{
              cds.layer <- NULL
            }

            # ==================================================================
            # PLOT
            # whether add new signal range
            if(new_signal_range == TRUE){
              axis.text.y = element_blank()
              axis.ticks.y = element_blank()

              range <- paste("[0-",round(max(pldf$rpm),digits = range_digit),"]",sep = "")

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

            # whether show merged sample and seperate sample together
            if(sep_mer_sample == TRUE){
              tmp.df2 <- pldf
              tmp.df2$sp <- "merged sample"

              pldf$sp <- pldf$sample
              pldf <- rbind(tmp.df2, pldf)

              facet_layer <- ggplot2::facet_grid(sp~rname,switch = "y")
            }

            # sample orders
            if(!is.null(sample_order)){
              if(sep_mer_sample == TRUE){
                pldf$sp <- factor(pldf$sp ,levels = c(sample_order,"merged sample"))
              }else{
                pldf$sample <- factor(pldf$sample ,levels = sample_order)
              }

            }

            # ==================================================================
            # plot
            # p <-
            ggplot(pldf) +
              player +
              range_label +
              theme_bw() +
              facet_layer +
              theme(panel.grid = element_blank(),
                    axis.text = element_text(colour = "black"),
                    strip.text = element_text(colour = "black",size = rel(1)),
                    strip.text.x = element_text(face = "bold",colour = "black",size = rel(1)),
                    strip.background = element_blank(),
                    strip.placement = "outside",
                    strip.text.y.left = element_text(angle = 0,hjust = 1),
                    axis.text.y = axis.text.y,
                    axis.ticks.y = axis.ticks.y,
                    ggside.panel.scale.x = structure_panel_height,
                    ggside.panel.background = element_blank(),
                    ggside.panel.border = element_blank()) +
              # gene structure
              bg.layer + exon.layer + cds.layer +
              ggside::scale_xsidey_continuous(breaks = NULL,limits = sturc_limit) +
              xlab("Position along genome (nt)") +
              ylab(ylab) +
              ggside(collapse = "x") +
              col + col1

          }
)

