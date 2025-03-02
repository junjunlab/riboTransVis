# ==============================================================================
# function for length plot
# ==============================================================================

#' Generate Length Distribution Plot
#'
#' @description This generic function generates a length distribution plot for sequencing reads.
#'
#' @param object An object containing experimental data.
#' @param ... Additional arguments passed to class-specific methods.
#'
#' @return A ggplot object representing the length distribution of reads.
#'
#' @export
#' @rdname length_plot
setGeneric("length_plot",function(object,...) standardGeneric("length_plot"))



#' @title Length Distribution Plot for Ribo-seq Reads
#' @description The `length_plot` function visualizes the distribution of read lengths in a Ribo-seq dataset.
#' It can either show the total counts of different read lengths (`type = "length"`)
#' or display frame-periodicity information (`type = "frame_length"`).
#'
#' @param object A `ribotrans` object that contains summary read count information.
#' @param read_length A numeric vector of length 2 (default: `c(20, 35)`). Defines the range of read lengths to visualize.
#' @param text_size Numeric (default: `2.5`). Specifies the font size for periodicity labels, applicable when `type = "frame_length"`.
#' @param add_periodicity_label Logical (default: `TRUE`). If `TRUE`, adds periodicity percentage labels to the plot (only for `type = "frame_length"`).
#' @param return_data Logical (default: `FALSE`). If `TRUE`, returns the processed data frame instead of a plot.
#' @param type Character string, either `"length"` (default) or `"frame_length"`.
#'   - `"length"`: Displays a bar plot showing the distribution of read lengths.
#'   - `"frame_length"`: Displays a bar plot with frame-periodicity information, highlighting the reading frame (0, 1, 2).
#'
#' @return A `ggplot2` object showing the read length distribution or a processed data frame if `return_data = TRUE`.
#'
#' @details
#' - If `type = "length"`, the function summarizes read counts by length.
#' - If `type = "frame_length"`, the function calculates frame-periodicity and shows the proportion of reads mapped to each frame.
#' - Frame-periodicity is computed as: ***percentage of reads in frame-0 relative to all frames***.
#' - Reads outside annotated coding regions (`mstart !=0` or `mstop !=0`) are analyzed for periodicity.
#'
#' @examples
#' \dontrun{
#' data(ribo_obj)  # Assume ribo_obj is a 'ribotrans' object
#'
#' # Plot read length distribution
#' length_plot(ribo_obj, read_length = c(25, 30), type = "length")
#'
#' # Plot frame-periodicity with periodicity labels
#' length_plot(ribo_obj, type = "frame_length", add_periodicity_label = TRUE)
#'
#' # Retrieve processed data
#' df <- length_plot(ribo_obj, return_data = TRUE, type = "frame_length")
#' head(df)
#' }
#'
#' @seealso [ggplot2::geom_col()], [ggplot2::facet_wrap()], [dplyr::mutate()], [scales::label_log()]
#'
#'
#' @importFrom ggplot2 ggplot geom_col scale_y_continuous facet_wrap
#' @importFrom ggplot2 theme element_text element_blank position_dodge2
#' @importFrom scales label_log
#' @importFrom fastplyr f_group_by f_summarise f_filter f_select
#'
#' @export
#' @rdname length_plot
setMethod("length_plot",
          signature(object = "ribotrans"),
          function(object,
                   read_length = c(20,35),
                   text_size = 2.5,
                   add_periodicity_label = TRUE,
                   return_data = FALSE,
                   type = c("length","frame_length")){
            type <- match.arg(type, choices = c("length","frame_length"))

            # check plot type
            if(type == "length"){
              pltdf <- object@summary_info %>%
                fastplyr::f_group_by(sample, qwidth) %>%
                fastplyr::f_summarise(counts = sum(count)) %>%
                dplyr::filter(qwidth >= read_length[1] & qwidth <= read_length[2])

              ptlayer <- geom_col(aes(x = qwidth, y = counts),width = 0.9,fill = "grey30")
              periodicity.layer <- NULL
              col <- NULL
            }else{
              pltdf <- object@summary_info %>%
                fastplyr::f_filter(mstart != 0 | mstop != 0) %>%
                dplyr::mutate(frame = (pos - mstart)%%3) %>%
                fastplyr::f_select(sample,qwidth,frame,count) %>%
                fastplyr::f_group_by(sample, qwidth, frame) %>%
                fastplyr::f_summarise(counts = sum(count)) %>%
                dplyr::filter(qwidth >= read_length[1] & qwidth <= read_length[2])

              # calculate periodicy
              read.length.all <- pltdf %>%
                fastplyr::f_group_by(sample, qwidth) %>%
                fastplyr::f_summarise(all_counts = sum(counts))

              pltdf <- pltdf %>%
                dplyr::left_join(y = read.length.all,by = c("sample", "qwidth")) %>%
                dplyr::mutate(periodicity = round((counts/all_counts)*100,digits = 1))

              ptlayer <- geom_col(aes(x = qwidth, y = counts, fill = factor(frame)),
                                  width = 0.9, position = position_dodge2())

              if(add_periodicity_label == TRUE){
                periodicity.layer <- geom_text(data = pltdf %>% dplyr::filter(frame == 0),
                                               aes(x = qwidth,y = max(counts),label = periodicity),
                                               size = text_size)
              }else{
                periodicity.layer <- NULL
              }

              col <- scale_fill_manual(values = c("0" = "#003366", "1" = "#336699", "2" = "#CCCCCC"),
                                       name = "frame")
            }

            # plot
            p <-
              ggplot(pltdf) +
              ptlayer +
              periodicity.layer +
              col +
              theme(panel.grid = element_blank(),
                    strip.text = element_text(colour = "black",face = "bold",size = rel(1)),
                    axis.text = element_text(colour = "black")) +
              facet_wrap(~sample,scales = "free") +
              scale_y_continuous(labels = scales::label_log(base = 10,digits = 1)) +
              xlab("Read length (nt)") + ylab("Number of reads")

            # return
            if(return_data == FALSE){
              return(p)
            }else{
              return(pltdf)
            }
          }
)




# ==============================================================================
# function for relative distance plot
# ==============================================================================

#' Generate Relative Distance Plot
#'
#' @description This generic function generates a relative distance distribution plot for sequencing reads.
#'
#' @param object An object containing experimental data.
#' @param ... Additional arguments passed to class-specific methods.
#'
#' @return A `ggplot2` object representing the relative distance distribution of reads.
#'
#' @export
setGeneric("relative_dist_plot",function(object,...) standardGeneric("relative_dist_plot"))



#' @describeIn relative_dist_plot Generate a relative distance plot for `ribotrans` objects.
#'
#' @description This method generates a plot showing the relative distance of reads to either
#' the start codon or stop codon, grouped by frame.
#'
#' @param object A `ribotrans` object containing summary information.
#' @param return_data Whether return the data for plot, default `FALSE`.
#' @param type Character string specifying the plot type. Must be one of:
#' \itemize{
#'   \item `"rel2start"` - Relative to the start codon.
#'   \item `"rel2stop"` - Relative to the stop codon.
#' }
#' @param read_length Numeric vector of length two specifying the filter range for read length (default: `c(20,35)`).
#' @param rel2st_dist Numeric vector of length two specifying the relative distance range to the start codon (default: `c(-50,100)`).
#' @param rel2sp_dist Numeric vector of length two specifying the relative distance range to the stop codon (default: `c(-100,50)`).
#' @param facet_wrap A `ggh4x::facet_wrap2()` object for facetting by sample (default behavior).
#' @param facet_grid A `ggh4x::facet_grid2()` object for facetting by read length and sample.
#' @param collapse_readlength Logical, whether to collapse the visualization across all read lengths (`TRUE`) or show separate read lengths (`FALSE`, default).
#'
#' @return A `ggplot2` object visualizing the relative distribution to start/stop codons by frame.
#'
#' @details
#' - If `type = "rel2start"`, the distance is calculated as `pos - mstart` (relative to start codon).
#' - If `type = "rel2stop"`, the distance is calculated as `pos - mstop` (relative to stop codon).
#' - Filtering is applied based on `read_length`, `rel2st_dist`, and `rel2sp_dist`.
#' - The plot can either group data by read length (`collapse_readlength = FALSE`) or aggregate reads (`collapse_readlength = TRUE`).
#'
#' @examples
#' \dontrun{
#' # Generate relative plot with start codon distances
#' relative_dist_plot(ribo_obj, type = "rel2start")
#'
#' # Generate relative plot with stop codon distances
#' relative_dist_plot(ribo_obj, type = "rel2stop")
#' }
#'
#' @importFrom ggplot2 ggplot geom_col scale_fill_manual facet_wrap
#' @importFrom ggplot2 theme element_text element_blank position_dodge2
#' @importFrom scales label_log
#' @importFrom ggh4x facet_wrap2 facet_grid2
#' @importFrom fastplyr f_group_by f_summarise f_filter f_select
#' @importFrom dplyr mutate filter
#'
#' @export
setMethod("relative_dist_plot",
          signature(object = "ribotrans"),
          function(object,
                   return_data = FALSE,
                   type = c("rel2start","rel2stop"),
                   read_length = c(20,35),
                   rel2st_dist = c(-50,100),
                   rel2sp_dist = c(-100,50),
                   facet_wrap = ggh4x::facet_wrap2(~sample,scales = "free"),
                   facet_grid = ggh4x::facet_grid2(qwidth~sample,scales = "free", independent = "y"),
                   collapse_readlength = FALSE){
            type <- match.arg(type, choices = c("rel2start","rel2stop"))

            # extact data
            # check plot type
            if(type == "rel2start"){
              pltdf <- object@summary_info %>%
                fastplyr::f_filter(mstart != 0 | mstop != 0) %>%
                dplyr::mutate(frame = (pos - mstart)%%3,
                              rel = pos - mstart) %>%
                fastplyr::f_select(sample,qwidth,frame,rel,count)
            }else{
              pltdf <- object@summary_info %>%
                fastplyr::f_filter(mstart != 0 | mstop != 0) %>%
                dplyr::mutate(frame = (pos - mstart)%%3,
                              rel = pos - mstop) %>%
                fastplyr::f_select(sample,qwidth,frame,rel,count)
            }

            # filter length and distance
            if(type == "rel2start"){
              pltdf <- pltdf %>%
                dplyr::filter(qwidth >= read_length[1] & qwidth <= read_length[2]) %>%
                dplyr::filter(rel >= rel2st_dist[1] & rel <= rel2st_dist[2])
            }else{
              pltdf <- pltdf %>%
                dplyr::filter(qwidth >= read_length[1] & qwidth <= read_length[2]) %>%
                dplyr::filter(rel >= rel2sp_dist[1] & rel <= rel2sp_dist[2])
            }

            # check whether plot for each read length
            if(collapse_readlength == TRUE){
              pltdf <- pltdf %>%
                fastplyr::f_group_by(sample, rel, frame) %>%
                fastplyr::f_summarise(counts = sum(count))

              facetlayer <- facet_wrap
            }else{
              pltdf <- pltdf %>%
                fastplyr::f_group_by(sample,qwidth, rel, frame) %>%
                fastplyr::f_summarise(counts = sum(count))

              facetlayer <- facet_grid
            }


            # plot
            p <-
              ggplot(pltdf) +
              geom_col(aes(x = rel, y = counts, fill = factor(frame)),
                       width = 0.9, position = position_dodge2()) +
              theme(panel.grid = element_blank(),
                    strip.text = element_text(colour = "black",face = "bold",size = rel(1)),
                    axis.text = element_text(colour = "black")) +
              # facet_grid(qwidth~sample,scales = "free") +
              facetlayer +
              scale_y_continuous(labels = scales::label_log(base = 10,digits = 1)) +
              scale_fill_manual(values = c("0" = "#003366", "1" = "#336699", "2" = "#CCCCCC"),
                                name = "frame") +
              xlab("Distance to start/stop codon (nt)") + ylab("Number of reads")

            # return
            if(return_data == FALSE){
              return(p)
            }else{
              return(pltdf)
            }
          }
)



# ==============================================================================
# function for relative heatmap plot
# ==============================================================================

#' Relative Heatmap Plot Generic Function
#'
#' @description This is a generic function for creating relative heatmap plots
#' across different data objects.
#'
#' @param object An object that contains the necessary data for plotting.
#' @param ... Additional arguments passed to specific methods.
#'
#' @return A heatmap plot specific to the method invoked.
#'
#' @export
#' @rdname relative_heatmap_plot
setGeneric("relative_heatmap_plot",function(object,...) standardGeneric("relative_heatmap_plot"))



#' Generate a Relative Heatmap Plot
#'
#' @description This function creates a heatmap that represents the relative position
#' of reads with respect to either start or stop codons of genes.
#'
#'
#' @param object An object of class \code{ribotrans}, containing summary information from ribosome profiling data.
#' @param type A character string indicating the type of plot. Options are:
#'  \code{"rel2start"} to plot relative to start codon or
#'  \code{"rel2stop"} to plot relative to stop codon. Default is \code{"rel2start"}.
#' @param read_length A numeric vector of length 2 specifying the range of read lengths to include in the plot. Default is \code{c(20, 35)}.
#' @param rel_dist A numeric vector of length 2 specifying the range of relative distances to plot. Default is \code{c(-100, 100)}.
#' @param log_scale A logical value indicating whether to apply a logarithmic transformation to counts. Default is \code{FALSE}.
#' @param return_data A logical value determining the return type. If \code{TRUE}, the summarized data is returned; if \code{FALSE}, the plot object is returned. Default is \code{FALSE}.
#'
#' @return A ggplot object if \code{return_data} is \code{FALSE}, otherwise a data.table with summarized counts and relevant metadata.
#'
#' @details
#' The function checks whether the summary information is available in the provided
#' \code{ribotrans} object and filters data based on the specified criteria.
#' The resulting heatmap visualizes the distribution of read counts along with
#' their relative positions to the coding features.
#'
#' @import ggplot2
#' @importFrom dplyr mutate filter group_by summarise
#' @importFrom fastplyr f_filter f_group_by f_summarise
#'
#'
#' @export
#' @rdname relative_heatmap_plot
setMethod("relative_heatmap_plot",
          signature(object = "ribotrans"),
          function(object,
                   type = c("rel2start","rel2stop"),
                   read_length = c(20,35),
                   rel_dist = c(-100,100),
                   log_scale = FALSE,
                   return_data = FALSE){

            type <- match.arg(type,choices = c("rel2start","rel2stop"))

            # check data
            if(nrow(object@summary_info) == 0){
              stop("Please run `generate_summary` first!")
            }

            # extarct data to plot
            if(type == "rel2start"){
              summary.info <- object@summary_info %>%
                fastplyr::f_filter(mstart > 0 & mstop > 0) %>%
                dplyr::mutate(rel_pos = pos - mstart)

              xlab <- "Distance to start codon (nt)"
            }else{
              summary.info <- object@summary_info %>%
                fastplyr::f_filter(mstart > 0 & mstop > 0) %>%
                dplyr::mutate(rel_pos = pos - mstop)

              xlab <- "Distance to stop codon (nt)"
            }

            # filter data
            summary.info <- summary.info %>%
              fastplyr::f_filter(rel_pos >= rel_dist[1] & rel_pos <= rel_dist[2]) %>%
              fastplyr::f_filter(qwidth >= read_length[1] & qwidth <= read_length[2]) %>%
              fastplyr::f_group_by(sample, rel_pos, qwidth) %>%
              fastplyr::f_summarise(counts = sum(count))

            # check log_scale
            if(log_scale == TRUE){
              layer <- geom_tile(aes(x = rel_pos,y = qwidth, fill = log2(counts)))
            }else{
              layer <- geom_tile(aes(x = rel_pos,y = qwidth, fill = counts))
            }

            # plot
            p <-
              ggplot(summary.info) +
              layer +
              facet_wrap(~sample,scales = "free") +
              theme(axis.text = element_text(colour = "black"),
                    panel.grid = element_blank(),
                    strip.text = element_text(colour = "black",size = rel(1),face = "bold")) +
              scale_fill_viridis_c(option = "mako") +
              xlab(xlab) +
              ylab("Reads length (nt)")


            # check return type
            if(return_data == FALSE){
              return(p)
            }else{
              return(summary.info)
            }

          }
)




# ==============================================================================
# function for reads offset check plot
# ==============================================================================
#' @title Generic function for plotting relative offset in ribosome profiling data
#'
#' @description This generic function is used to plot the relative offset of reads
#' in ribosome profiling data.
#'
#' @param object An object of class \code{ribotrans}.
#' @param ... Additional arguments passed to specific methods.
#'
#' @return A plot showing the relative offset of reads.
#'
#' @seealso \code{\link{relative_offset_plot,ribotrans-method}}
#'
#' @export
#' @rdname relative_offset_plot
setGeneric("relative_offset_plot",function(object,...) standardGeneric("relative_offset_plot"))



#' @title Plot relative offset distribution in ribosome profiling data
#'
#' @description This method generates a plot of relative offsets of ribosome profiling
#' reads mapped to either the start or stop codon.
#'
#' @param object A \code{ribotrans} object containing ribosome profiling data.
#' @param type Character string specifying the reference point: \code{"rel2start"} for
#' start codon or \code{"rel2stop"} for stop codon. Default is \code{"rel2start"}.
#' @param read_length A numeric vector of length two specifying the range of read
#' lengths to include (default: \code{c(20, 35)}).
#' @param rel_dist A numeric vector of length two specifying the relative distance
#' range to include (default: \code{c(-30, 30)}).
#' @param rm_yaxis_label Logical; if \code{TRUE}, removes y-axis labels (default: \code{TRUE}).
#' @param max_offset_labelx Numeric value (between 0 and 1) specifying the x position
#' of the label showing the most frequent offset (default: \code{0.9}).
#' @param max_offset_labely Numeric value (between 0 and 1) specifying the y position
#' of the label showing the most frequent offset (default: \code{0.9}).
#' @param return_offset Logical; if \code{TRUE}, returns the most frequently observed
#' relative offsets instead of plotting (default: \code{FALSE}).
#'
#' @details This function visualizes the distribution of ribosome footprint
#' offsets relative to the start or stop codon. The plot includes read length on
#' the y-axis and relative position on the x-axis. The most frequent offset positions
#' are labeled, and vertical dashed lines highlight expected P-site positions.
#'
#' @return If \code{return_offset = FALSE}, returns a \code{ggplot} object with
#' the relative offset plot. If \code{return_offset = TRUE}, returns a data frame
#' containing the most frequently observed offsets.
#'
#' @examples
#' \dontrun{
#'   # Generate summary information before plotting
#'   ribo_obj <- generate_summary(ribo_obj)
#'
#'   # Plot relative offsets with default parameters
#'   relative_offset_plot(ribo_obj)
#'
#'   # Get the most frequent offset information
#'   offsets <- relative_offset_plot(ribo_obj, return_offset = TRUE)
#' }
#'
#' @seealso \code{\link{generate_summary}}
#'
#' @export
#' @rdname relative_offset_plot
setMethod("relative_offset_plot",
          signature(object = "ribotrans"),
          function(object,
                   type = c("rel2start","rel2stop"),
                   read_length = c(20,35),
                   rel_dist = c(-30,30),
                   rm_yaxis_label = TRUE,
                   max_offset_labelx = 0.9,
                   max_offset_labely = 0.9,
                   return_offset = FALSE){

            type <- match.arg(type,choices = c("rel2start","rel2stop"))

            # check data
            if(nrow(object@summary_info) == 0){
              stop("Please run `generate_summary` first!")
            }

            # extarct data to plot
            if(type == "rel2start"){
              summary.info <- object@summary_info %>%
                fastplyr::f_filter(mstart > 0 & mstop > 0) %>%
                dplyr::mutate(rel_pos = pos - mstart)

              xlab <- "Distance to start codon (nt)"
              xpos <- c(-12,-15,-18)
            }else{
              summary.info <- object@summary_info %>%
                fastplyr::f_filter(mstart > 0 & mstop > 0) %>%
                dplyr::mutate(rel_pos = pos - mstop)

              xlab <- "Distance to stop codon (nt)"
              xpos <- c(12,15,18)
            }

            # filter data
            summary.info <- summary.info %>%
              fastplyr::f_filter(rel_pos >= rel_dist[1] & rel_pos <= rel_dist[2]) %>%
              fastplyr::f_filter(qwidth >= read_length[1] & qwidth <= read_length[2]) %>%
              fastplyr::f_group_by(sample, rel_pos, qwidth) %>%
              fastplyr::f_summarise(counts = sum(count))

            # check y axis
            if(rm_yaxis_label == TRUE){
              axis.text.y = element_blank()
              axis.ticks.y = element_blank()
            }else{
              axis.text.y = NULL
              axis.ticks.y = NULL
            }

            # max offset position
            mx.pos <- summary.info %>%
              dplyr::slice_max(order_by = counts,n = 1,by = c(sample,qwidth))

            if (requireNamespace("ggpp", quietly = TRUE)) {
              mx.label <- ggpp::geom_label_npc(data = mx.pos,aes(npcx = max_offset_labelx,
                                                                 npcy = max_offset_labely,
                                                                 label = rel_pos))
            } else {
              warning("Package 'ggpp' is needed for this function to work.")
              mx.label <- NULL
            }


            # plot
            p <-
              ggplot(summary.info) +
              geom_path(aes(x = rel_pos,y = counts)) +
              mx.label +
              geom_vline(xintercept = xpos,lty ="dashed", color = "orange") +
              ggh4x::facet_grid2(qwidth~sample,scales = "free",independent = "y",axes = "y") +
              theme(axis.text = element_text(colour = "black"),
                    panel.grid = element_blank(),
                    axis.text.y = axis.text.y,
                    axis.ticks.y = axis.ticks.y,
                    strip.text = element_text(colour = "black",size = rel(1),face = "bold")) +
              scale_fill_viridis_c(option = "mako") +
              scale_y_continuous(labels = scales::label_log(base = 10,digits = 1)) +
              xlab(xlab) +
              ylab("Reads length (nt)")


            # check return type
            if(return_offset == FALSE){
              return(p)
            }else{
              return(mx.pos)
            }

          }
)

# ==============================================================================
# function for metagene plot
# ==============================================================================

#' Generate Metagene Profile Plot
#'
#' @description A generic function to generate a metagene profile plot
#' visualizing read distributions relative to start or stop codons.
#'
#' @param object An object containing sequencing data.
#' @param ... Additional arguments passed to specific methods.
#'
#' @return A `ggplot2` object showing the metagene profile distribution.
#'
#' @export
setGeneric("metagene_plot",function(object,...) standardGeneric("metagene_plot"))



#' @describeIn metagene_plot Generate a metagene profile plot for `ribotrans` objects.
#'
#' @description This method creates a metagene profile plot showing the
#' distribution of ribosome-protected fragments (RPFs) relative to the
#' start or stop codon.
#'
#' @param object A `ribotrans` object containing RNA sequencing metrics.
#' @param type Character string specifying the metagene type. Must be one of:
#' \itemize{
#'   \item `"rel2start"` - Distance relative to the start codon.
#'   \item `"rel2stop"` - Distance relative to the stop codon.
#' }
#' @param smooth Logical indicating if a rolling average should be applied to smooth the data (default: `FALSE`).
#' @param return_data Logical, if `TRUE`, returns the processed data instead of a plot (default: `FALSE`).
#' @param show_frame Logical, if `TRUE`, visualizes data by reading frame (default: `FALSE`).
#' @param read_length Numeric vector of length two specifying the range of read lengths to include (default: `c(25,31)`).
#' @param slide_window Integer defining the sliding window size for smoothing when `smooth = TRUE` (default: `3`).
#' @param rel2st_dist Numeric vector of length two, specifying distance range to include relative to start codon (default: `c(-50,100)`).
#' @param rel2sp_dist Numeric vector of length two, specifying distance range to include relative to stop codon (default: `c(-100,50)`).
#' @param facet_wrap A `ggplot2::facet_wrap()` object to control faceting by sample (default behavior).
#'
#' @return If `return_data = FALSE`, returns a `ggplot2` object showing the metagene profile plot.
#' If `return_data = TRUE`, returns a processed `data.frame` containing metagene read distributions.
#'
#' @details
#' - The function normalizes read counts based on total mapped reads (RPM normalization).
#' - Optionally, it applies a rolling average (`slide_window`) to smooth the profile.
#' - If `show_frame = TRUE`, colors are applied for individual reading frames.
#'
#' @examples
#' \dontrun{
#' # Generate standard metagene plot relative to start codon
#' metagene_plot(ribo_obj, type = "rel2start")
#'
#' # Generate smoothed metagene plot
#' metagene_plot(ribo_obj, type = "rel2stop", smooth = TRUE, slide_window = 3)
#'
#' # Get the underlying data instead of a plot
#' data <- metagene_plot(ribo_obj, type = "rel2start", return_data = TRUE)
#' head(data)
#' }
#'
#' @importFrom ggplot2 ggplot geom_line scale_color_manual facet_wrap
#' @importFrom ggplot2 theme element_text element_blank
#' @importFrom dplyr mutate filter left_join select summarise
#' @importFrom fastplyr f_group_by f_summarise f_filter f_select f_left_join
#' @importFrom purrr map_df
#'
#' @export
setMethod("metagene_plot",
          signature(object = "ribotrans"),
          function(object,
                   type = c("rel2start","rel2stop"),
                   smooth = FALSE,
                   return_data = FALSE,
                   show_frame = FALSE,
                   read_length = c(25,31),
                   slide_window = 3,
                   rel2st_dist = c(-50,100),
                   rel2sp_dist = c(-100,50),
                   facet_wrap = ggplot2::facet_wrap(~sample)){
            type <- match.arg(type, choices = c("rel2start","rel2stop"))

            features <- object@features

            pltdf <- object@summary_info %>%
              fastplyr::f_filter(mstart != 0 | mstop != 0,
                                 qwidth >= read_length[1] & qwidth <= read_length[2]) %>%
              fastplyr::f_group_by(sample,rname, pos) %>%
              fastplyr::f_summarise(counts = sum(count))

            # get total mapped reads
            lib <- object@library
            dpt <- subset(lib, type == "ribo" & sample %in% unique(pltdf$sample)) %>%
              dplyr::select(mappped_reads,sample)

            # rpm normalization
            pltdf <- pltdf %>%
              dplyr::left_join(y = dpt,by = "sample") %>%
              dplyr::mutate(rpm = (counts/mappped_reads)*10^6) %>%
              dplyr::select(-mappped_reads)

            # total transcipts
            ttt <- length(unique(pltdf$rname))

            pltdf <- pltdf %>%
              fastplyr::f_left_join(y = features[,c("idnew","mstart","mstop")],
                                    by = c("rname" = "idnew"))


            # calculate frame and relative distance
            if(type == "rel2start"){
              pltdf <- pltdf %>%
                dplyr::mutate(frame = (pos - mstart)%%3,
                              rel = pos - mstart) %>%
                fastplyr::f_select(sample,frame,pos,rel,rpm)
            }else{
              pltdf <- pltdf %>%
                dplyr::mutate(frame = (pos - mstop)%%3,
                              rel = pos - mstop) %>%
                fastplyr::f_select(sample,frame,pos,rel,rpm)
            }


            # check whether plot for each read length
            if(show_frame == TRUE){
              pltdf <- pltdf %>%
                fastplyr::f_group_by(sample,rel,frame) %>%
                fastplyr::f_summarise(normval = sum(rpm)/ttt)

              pltlayer <- geom_line(aes(x = rel, y = relexp, color = factor(frame)))
              col <- scale_fill_manual(values = c("0" = "#003366", "1" = "#336699", "2" = "#CCCCCC"),
                                       name = "frame")
            }else{
              pltdf <- pltdf %>%
                fastplyr::f_group_by(sample,rel) %>%
                fastplyr::f_summarise(normval = sum(rpm)/ttt)

              pltlayer <- geom_line(aes(x = rel, y = relexp))
              col <- NULL
            }


            # filter length and distance
            if(type == "rel2start"){
              pltdf <- pltdf %>%
                dplyr::filter(rel >= rel2st_dist[1] & rel <= rel2st_dist[2])

              distrg <- seq(rel2st_dist[1],rel2st_dist[2],1)
            }else{
              pltdf <- pltdf %>%
                dplyr::filter(rel >= rel2sp_dist[1] & rel <= rel2sp_dist[2])

              distrg <- seq(rel2sp_dist[1],rel2sp_dist[2],1)
            }

            # ==================================================================
            # whether smooth data
            if(smooth == TRUE){
              sp <- unique(pltdf$sample)

              # x = 1
              purrr::map_df(seq_along(sp),function(x){
                tmp <- subset(pltdf, sample == sp[x])

                tmp2 <- data.frame(sample = sp[x],rel = distrg)
                tmp2 <- tmp2 %>%
                  dplyr::left_join(y = tmp,by = c("sample","rel"))

                if (requireNamespace("zoo", quietly = TRUE)) {
                  tmp2$smooth <- zoo::rollmean(tmp2$normval, k = slide_window, fill = 0)
                } else {
                  warning("Package 'zoo' is needed for this function to work.")
                }

                return(tmp2)
              })-> sm.df
            }else{
              sm.df <- pltdf %>%
                dplyr::mutate(smooth = normval)

            }


            ave.tmp <- sm.df %>%
              dplyr::group_by(sample) %>%
              dplyr::summarise(allexp = sum(smooth)) %>%
              dplyr::mutate(ave = allexp/length(distrg)) %>%
              dplyr::select(sample, ave)

            sm.df <- sm.df %>%
              dplyr::left_join(y = ave.tmp,by = "sample") %>%
              dplyr::mutate(relexp = smooth/ave)

            # plot
            p <-
              ggplot(sm.df) +
              pltlayer +
              theme(panel.grid = element_blank(),
                    strip.text = element_text(colour = "black",face = "bold",size = rel(1)),
                    axis.text = element_text(colour = "black")) +
              # facet_wrap(~sample) +
              facet_wrap +
              col +
              xlab("Distance to start/stop codon (nt)") + ylab("Normalized reads")


            # return
            if(return_data == FALSE){
              return(p)
            }else{
              return(sm.df)
            }

          }
)
