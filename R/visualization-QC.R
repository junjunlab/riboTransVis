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
# function for feature plot
# ==============================================================================

#' Plot Read Distribution Across RNA Features
#'
#' This function generates a bar plot showing the distribution of ribosome profiling reads
#' across different RNA features (5' UTR, CDS, 3' UTR).
#'
#' @param object A `ribotrans` object that contains summary read information.
#' @param return_data A logical value indicating whether to return the processed
#' summary data instead of the plot. Default is `FALSE` (returns a `ggplot` object).
#' @param ... Additional arguments (currently unused).
#'
#' @return If `return_data = FALSE`, a `ggplot2` object visualizing the **read counts**
#' in different regions (_5' UTR, CDS, 3' UTR_) is returned.
#' If `return_data = TRUE`, a `data.frame` summarizing the read counts by region
#' and sample is returned.
#'
#' @details
#' - Reads are assigned to **5' UTR, CDS, or 3' UTR** based on the midpoint of the read position.
#' - The function uses `dplyr` and `fastplyr` for efficient data manipulation.
#' - The resulting bar plot groups reads by `sample` and scales axes using a log-10 transformation.
#'
#' @examples
#' \dontrun{
#' # Example: Generate a feature plot
#' feature_plot(obj)
#'
#' # Example: Retrieve summarized data
#' res <- feature_plot(obj, return_data = TRUE)
#' head(res)
#' }
#'
#' @export
setGeneric("feature_plot",function(object,...) standardGeneric("feature_plot"))


#' @rdname feature_plot
#' @export
setMethod("feature_plot",
          signature(object = "ribotrans"),
          function(object,return_data = FALSE){
            # assign features
            sry <- object@summary_info %>%
              fastplyr::f_filter(mstart != 0 | mstop != 0) %>%
              dplyr::mutate(type = dplyr::case_when(pos + qwidth/2 < mstart ~ "5'UTR",
                                                    pos + qwidth/2 >= mstart & pos + qwidth/2 < mstop ~ "CDS",
                                                    pos + qwidth/2 >= mstop ~ "3'UTR")) %>%
              fastplyr::f_group_by(sample, type) %>%
              fastplyr::f_summarise(counts = sum(count))

            # plot
            # order
            sry$type <- factor(sry$type,levels = c("CDS","5'UTR","3'UTR"))

            p <-
              ggplot(sry) +
              geom_col(aes(x = type,y = counts,fill = type),width = 0.6,show.legend = F) +
              facet_wrap(~sample,scales = "free_y")+
              theme(panel.grid = element_blank(),
                    strip.text = element_text(colour = "black",face = "bold",size = rel(1)),
                    axis.text = element_text(colour = "black")) +
              scale_y_continuous(labels = scales::label_log(base = 10,digits = 1)) +
              xlab("Read features") + ylab("Number of reads")

            # return
            if(return_data == FALSE){
              return(p)
            }else{
              return(sry)
            }
          }
)



# ==============================================================================
# function for whole_metagene plot
# ==============================================================================

#' Generate Metagene Profile Plot for Ribosome Density
#'
#' This function creates a **metagene profile plot** to visualize ribosome density
#' across different transcript features (5' UTR, CDS, and 3' UTR).
#'
#' @param object A `ribotrans` object containing summarized ribosome profiling data.
#' @param auto_scale Logical. If `TRUE`, the function automatically scales UTR lengths
#' relative to the CDS region using the median UTR/CDS lengths across all transcripts.
#' Default is `FALSE`, meaning manual scaling is used (via `custom_scale`).
#' @param custom_scale A numeric vector of length **2** specifying the relative scaling
#' factors for **5' UTR and 3' UTR** (relative to CDS length).
#' **Default**: `c(0.3, 0.3)`, where both UTRs are **30% the size of the CDS region**.
#' @param geom_density_bw Bandwidth for `geom_density()`. Controls the smoothness
#' of the density plot. **Default**: `0.0005`.
#' @param geom_density_n Number of points used to generate the density estimate.
#' **Default**: `512`. Lower values can speed up computations.
#' @param return_data Logical. If `FALSE`, returns a `ggplot2` object for visualization.
#' If `TRUE`, returns **processed transcript-level ribosome density data** as a `data.frame`.
#' **Default**: `FALSE` (returns plot).
#' @param ... Additional arguments (currently unused).
#'
#' @return If `return_data = FALSE`, returns a `ggplot2` object displaying **ribosome density profiles** across `5' UTR`, `CDS`, and `3' UTR`.
#' If `return_data = TRUE`, returns a `data.frame` containing the **scaled relative ribosome density positions**.
#'
#' @details
#' - This function calculates a **relative distance metric** for each read's mapped position and scales them appropriately into **5' UTR, CDS, and 3' UTR regions**.
#' - If `auto_scale = TRUE`, the function automatically estimates median UTR/CDS
#' scaling factors from the dataset.
#' - The **X-axis** represents normalized transcript positions:
#'   - `1` → start codon of the **CDS**
#'   - `2` → stop codon of the **CDS**
#'   - `Scaled values between 0 and 1` → **5' UTR**
#'   - `Scaled values between 2 and 3` → **3' UTR**
#' - **Dashed reference lines** (`geom_vline()`) mark **CDS start (1)** and **CDS stop (2)**.
#'
#' @examples
#' \dontrun{
#' # Generate a metagene plot
#' whole_metagene_plot(obj)
#'
#' # Return processed density data instead of plot
#' density_data <- whole_metagene_plot(obj, return_data = TRUE)
#' head(density_data)
#'
#' # Adjust kernel density bandwidth for smoother plot
#' whole_metagene_plot(obj, geom_density_bw = 0.001)
#' }
#'
#' @export
setGeneric("whole_metagene_plot",function(object,...) standardGeneric("whole_metagene_plot"))




#' @rdname whole_metagene_plot
#' @export
setMethod("whole_metagene_plot",
          signature(object = "ribotrans"),
          function(object,
                   auto_scale = FALSE,
                   custom_scale = c(0.3,0.3),
                   geom_density_bw = 0.0005,
                   geom_density_n = 512,
                   return_data = FALSE){
            # ==================================================================
            # calculate relative distance
            sry <- object@summary_info %>%
              fastplyr::f_filter(mstart != 0 | mstop != 0) %>%
              dplyr::mutate(pos = pos + qwidth/2) %>%
              dplyr::mutate(rel = dplyr::case_when(pos < mstart ~ pos/mstart,
                                                   pos >= mstart & pos < mstop ~ (pos - mstart)/(mstop - mstart) + 1,
                                                   pos >= mstop ~ (pos - mstop)/(translen - mstop) + 2)) %>%
              fastplyr::f_select(sample, rel)

            # ==================================================================
            # scale to cds region
            fea <- object@features

            if(auto_scale == TRUE){
              if (requireNamespace("stats", quietly = TRUE)) {
                utr5 <- stats::median(fea$utr5)
                cds <- stats::median(fea$cds)
                utr3 <- stats::median(fea$utr3)
              } else {
                warning("Package 'stats' is needed for this function to work.")
              }

              sutr5 <- utr5/cds
              sutr3 <- utr3/cds
            }else{
              sutr5 <- custom_scale[1]
              sutr3 <- custom_scale[2]
            }

            # re-scale regions
            sry <- sry %>%
              dplyr::mutate(rel_scale = dplyr::case_when(
                rel < 1 ~ scales::rescale(rel,to = c(1 - sutr5,1),from = c(0,1)),
                rel >= 2 ~ scales::rescale(rel,to = c(2,2 + sutr3),from = c(2,3)),
                .default = rel))

            # ==================================================================
            # plot
            p <-
              ggplot(sry) +
              geom_density(aes(x = rel_scale),color = "#003366",bw = geom_density_bw, n = geom_density_n) +
              facet_wrap(~sample) +
              theme_bw() +
              theme(panel.grid = element_blank(),
                    strip.text = element_text(colour = "black",face = "bold",size = rel(1)),
                    axis.text = element_text(colour = "black")) +
              scale_x_continuous(breaks = c(1,1.5,2),
                                 limits = c(1 - sutr5,2 + sutr3),
                                 labels = c("start","CDS","stop")) +
              geom_vline(xintercept = c(1,2),lty = "dashed", color = "grey30") +
              xlab("Ribosome density (AU)")

            # return
            if(return_data == FALSE){
              return(p)
            }else{
              return(sry)
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
              offset <- mx.pos[,1:3]

              return(offset)
            }

          }
)

# ==============================================================================
# function for metagene plot
# ==============================================================================

#' Generate Metagene Plots for Ribosome Profiling Data
#'
#' This function creates **metagene plots**, visualizing **ribosome occupancy** relative
#' to **start/stop codons** in translated mRNA regions (CDS).
#'
#' @param object A `ribotrans` object containing ribosome profiling data.
#' @param do_offset_correct Logical. If `TRUE`, performs **offset correction**
#' using `do_offset_correction()`. **Default**: `FALSE`.
#' @param position_shift Integer defining how much to adjust **ribosome footprint positions**
#' during offset correction. **Default**: `0`.
#' @param type Character. Defines metagene analysis type:
#'   - `"rel2start"` (default): Plots occupancy relative to **start codon**.
#'   - `"rel2stop"`: Plots occupancy relative to **stop codon**.
#' @param return_data Logical. If `TRUE`, returns processed **occupancy data**
#' instead of plotting. **Default**: `FALSE`.
#' @param mode Character. Defines the **x-axis resolution**:
#'   - `"nt"` (default): Positions in **nucleotides (nt)**.
#'   - `"codon"`: Positions in **codon units (3-nt bins)**.
#' @param read_length Numeric vector, specifying accepted **ribosome-protected fragment (RPF) lengths**.
#' **Default**: `c(25,31)`.
#' @param rel2st_dist Numeric vector of length 2, defining the **window range** (**nt**)
#' for `"rel2start"` plots. **Default**: `c(-50,100)`.
#' @param rel2sp_dist Numeric vector of length 2, defining the **window range** (**nt**)
#' for `"rel2stop"` plots. **Default**: `c(-100,50)`.
#' @param facet_wrap A `ggplot2::facet_wrap()` object to define **faceting** for multiple samples.
#' **Default**: `ggplot2::facet_wrap(~sample)`.
#' @param ... Additional arguments (currently unused).
#'
#' @return
#' - If `return_data = FALSE`: Returns a `ggplot` object of the **metagene plot**.
#' - If `return_data = TRUE`: Returns a `data.frame` with:
#'   - `sample`: Sample names.
#'   - `rel`: Relative position (nt or codon).
#'   - `avg`: **Mean normalized ribosome occupancy** at each position.
#'
#' @details
#' - The function computes **ribosome footprint density** aligned to **start or stop codons**.
#' - If `mode = "codon"`, positions are converted to **codon units** (i.e., `rel/3`).
#' - The occupancy is **normalized by average translation density** to account for
#' differences in transcript lengths.
#'
#' @examples
#' \dontrun{
#' # Metagene profile aligned to start codon (default)
#' p <- metagene_plot(obj, type = "rel2start", mode = "nt")
#' print(p)
#'
#' # Return occupancy data instead of a plot
#' data <- metagene_plot(obj, type = "rel2stop", mode = "codon", return_data = TRUE)
#' head(data)
#'
#' # Customize plot with ggplot2
#' ggplot(data, aes(x = rel, y = avg, color = sample)) +
#'     geom_line() +
#'     theme_minimal()
#' }
#'
#' @importFrom Biostrings readDNAStringSet translate vmatchPattern
#' @importFrom ggplot2 ggplot geom_path theme element_blank element_text facet_wrap
#' @importFrom dplyr filter mutate select left_join rename
#' @importFrom fastplyr f_filter f_select f_group_by f_summarise
#' @importFrom purrr map_df
#' @export
setGeneric("metagene_plot",function(object,...) standardGeneric("metagene_plot"))



#' @rdname metagene_plot
#' @export
setMethod("metagene_plot",
          signature(object = "ribotrans"),
          function(object,
                   do_offset_correct = FALSE,
                   position_shift = 0,
                   type = c("rel2start","rel2stop"),
                   return_data = FALSE,
                   mode = c("nt", "codon"),
                   read_length = c(25,31),
                   rel2st_dist = c(-50,100),
                   rel2sp_dist = c(-100,50),
                   facet_wrap = ggplot2::facet_wrap(~sample)){
            type <- match.arg(type, choices = c("rel2start","rel2stop"))
            mode <- match.arg(mode, choices = c("nt", "codon"))

            # check type
            if(type == "rel2start"){
              dist <- rel2st_dist
              var <- "mstart"
            }else{
              dist <- rel2sp_dist
              var <- "mstop"
            }

            features <- object@features

            # whether do reads offset correction
            if(do_offset_correct == TRUE){
              sry <- do_offset_correction(object = object,shift = position_shift)
            }else{
              sry <- object@summary_info
            }

            sry <- sry %>% fastplyr::f_filter(mstart != 0 | mstop != 0)

            # average counts per position
            avg.ct <- sry %>%
              dplyr::mutate(cdslen = mstop - mstart + 1) %>%
              fastplyr::f_group_by(sample,rname,cdslen) %>%
              fastplyr::f_summarise(counts = sum(count)) %>%
              dplyr::mutate(avg_ct = counts/cdslen) %>%
              dplyr::select(sample,rname,avg_ct)

            # filter relative distance
            pltdf <- sry %>%
              dplyr::left_join(y = avg.ct,by = c("sample", "rname")) %>%
              dplyr::filter(mstop - mstart >= max(abs(dist))) %>%
              dplyr::mutate(rel = pos - .data[[var]], norm = count/avg_ct) %>%
              fastplyr::f_group_by(sample,rel) %>%
              fastplyr::f_summarise(normsm = sum(norm)) %>%
              dplyr::filter(rel >= dist[1] & rel <= dist[2])

            # average occupancy
            sm <- pltdf %>%
              fastplyr::f_group_by(sample) %>%
              fastplyr::f_summarise(norm_avg = sum(normsm)/(abs(dist[2]) + abs(dist[1]) + 1))

            pltdf <- pltdf %>%
              dplyr::left_join(y = sm,by = "sample") %>%
              dplyr::mutate(avg = normsm/norm_avg)

            # check show mode
            if(mode == "codon"){
              pltdf <- pltdf %>%
                dplyr::mutate(codon = rel %/% 3 + 1) %>%
                fastplyr::f_group_by(sample,codon) %>%
                fastplyr::f_summarise(avg = mean(avg)) %>%
                dplyr::rename(rel = codon)
            }

            # ==================================================================
            # plot
            p <-
              ggplot(pltdf) +
              geom_path(aes(x = rel,y = avg)) +
              theme(panel.grid = element_blank(),
                    strip.text = element_text(colour = "black",face = "bold",size = rel(1)),
                    axis.text = element_text(colour = "black")) +
              # facet_wrap(~sample) +
              facet_wrap +
              xlab("Distance to start/stop codon (nt)") +
              ylab("Average ribosome occupancy")

            # return
            if(return_data == FALSE){
              return(p)
            }else{
              return(pltdf)
            }
          }
)

