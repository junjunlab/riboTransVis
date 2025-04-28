# ==============================================================================
# get_occupancy for serp class
# ==============================================================================

#' Get Occupancy Data for SERP Objects
#'
#' @description
#' Extract and process occupancy data from ribosome profiling BAM files for an object of class 'serp'.
#' This method can process both total and IP (immunoprecipitation) samples.
#'
#' @param object An object of class 'serp'.
#' @param gene_name Character vector. Names of genes to extract data for. If NULL (default),
#'        data for all genes will be extracted.
#' @param smooth Logical. Whether to apply smoothing to the occupancy profile. Default: FALSE.
#' @param serp_exp Character. Experiment type to process, either "total" or "ip". Default: "total".
#' @param coordinate_to_trans Logical. Whether to convert genome coordinates to transcript coordinates
#'        when using genome mapping. Default: FALSE.
#' @param do_reads_offset Logical. Whether to apply read offset correction. Default: FALSE.
#' @param merge_rep Logical. If \code{TRUE}, average replicate samples within the same sample_group. Default: \code{FALSE}.
#' @param slide_window Numeric. Size of the sliding window for smoothing. Default: 30.
#' @param ... Additional arguments (currently unused).
#'
#' @return An updated 'serp' object with the occupancy data stored in the corresponding slot
#'         (either 'total_occupancy' or 'ip_occupancy').
#'
#' @details
#' This method extracts ribosome occupancy data from BAM files and processes it according
#' to the specified parameters. The method can handle both transcriptome and genome mappings,
#' and can apply various adjustments such as read offset correction and position smoothing.
#'
#' The processing workflow includes:
#' 1. Extracting occupancy data based on mapping type (genome or transcriptome)
#' 2. Optionally applying read offset correction based on stored offset information
#' 3. Aggregating counts and RPM (reads per million) values by position
#' 4. Optionally applying smoothing across positions
#' 5. Storing the processed data in the appropriate slot of the 'serp' object
#'
#' @examples
#' \dontrun{
#' # Create a serp object
#' serp_obj <- construct_serp(
#'   genome_file = "genome.fa",
#'   gtf_file = "annotations.gtf",
#'   total_bam_file = "total_sample.bam",
#'   total_sample_name = "total_rep1",
#'   IP_bam_file = "IP_sample.bam",
#'   IP_sample_name = "IP_rep1"
#' )
#'
#' # Get total occupancy data
#' serp_obj <- get_occupancy(serp_obj, serp_exp = "total")
#'
#' # Get IP occupancy data with smoothing
#' serp_obj <- get_occupancy(serp_obj,
#'                          serp_exp = "ip",
#'                          smooth = TRUE,
#'                          slide_window = 15)
#'
#' # Get occupancy for specific genes with read offset correction
#' serp_obj <- get_occupancy(serp_obj,
#'                          gene_name = "GENE1",
#'                          serp_exp = "total",
#'                          do_reads_offset = TRUE)
#' }
#'
#' @importFrom purrr map_df
#' @importFrom dplyr left_join mutate select group_by summarise
#' @importFrom stats na.omit
#'
#' @export
#' @rdname get_occupancy
setMethod("get_occupancy",
          signature(object = "serp"),
          function(object, gene_name = NULL,
                   smooth = FALSE,
                   serp_exp = c("total","ip"),
                   coordinate_to_trans = FALSE,
                   do_reads_offset = FALSE,
                   merge_rep = FALSE,
                   slide_window = 30){
            serp_exp <- match.arg(serp_exp,choices = c("total","ip"))

            assignment_mode <- object@assignment_mode

            bf <- subset(object@bam_file, type == serp_exp)
            ribobams <- bf$bam
            sp <- bf$sample
            gp <- bf$sample_group

            lib <- subset(object@library, type == serp_exp)
            lib <- lib[match(ribobams, lib$bam),]
            totalreads <- lib$mappped_reads

            # loop for each sample
            purrr::map_df(seq_along(ribobams),function(x){
              # check mapping type
              if(object@mapping_type == "genome"){
                tmp <- getOccupancyGenome(bam_file = ribobams[x],
                                          gene_name = gene_name,
                                          features = object@genome_trans_features,
                                          total_reads = totalreads[x],
                                          assignment_mode = assignment_mode,
                                          coordinate_to_trans = coordinate_to_trans)
              }else{
                tmp <- getOccupancy(bam_file = ribobams[x],
                                    gene_name = gene_name,
                                    features = object@features,
                                    total_reads = totalreads[x])
              }

              # add sample name
              tmp$sample <- sp[x]
              tmp$sample_group <- gp[x]

              # ================================================================
              # whether do read offset
              if(do_reads_offset == TRUE){
                # check offset data
                if(nrow(object@reads_offset_info) == 0){
                  stop("No reads offset information in object!")
                }

                offset <- object@reads_offset_info

                tmp$sample <- paste(sp[x], serp_exp, sep = "-")
                tmp$sample_group <- paste(gp[x], serp_exp, sep = "-")

                if(object@mapping_type == "genome" & coordinate_to_trans == FALSE){
                  adjusted <- tmp %>%
                    dplyr::left_join(y = offset[,c("sample","rel_pos","qwidth")],by = c("sample","qwidth")) %>%
                    stats::na.omit() %>%
                    dplyr::mutate(pos = dplyr::if_else(strand == "+",
                                                       pos - rel_pos,pos + rel_pos)) %>%
                    dplyr::select(-strand,-qwidth,-rel_pos)
                }else{
                  adjusted <- tmp %>%
                    dplyr::left_join(y = offset[,c("sample","rel_pos","qwidth")],by = c("sample","qwidth")) %>%
                    stats::na.omit() %>%
                    dplyr::mutate(pos = pos - rel_pos) %>%
                    dplyr::select(-qwidth,-rel_pos)
                }

                adjusted$sample <- sp[x]
                adjusted$sample_group <- gp[x]

              }else{
                adjusted <- tmp
              }

              adjusted <- adjusted %>%
                dplyr::group_by(sample,sample_group,rname,pos) %>%
                dplyr::summarise(count = sum(count),rpm = sum(rpm), .groups = "drop")

              return(adjusted)
            }) -> ribo.df

            # whether aggregate replicates
            if(merge_rep == TRUE){
              ribo.df <- ribo.df %>%
                fastplyr::f_group_by(sample_group,rname,pos) %>%
                fastplyr::f_summarise(count = ceiling(mean(count)),rpm = mean(rpm)) %>%
                dplyr::rename(sample = sample_group)
            }

            # smooth for each position
            cond <- object@mapping_type == "genome" & coordinate_to_trans == TRUE | object@mapping_type == "transcriptome"
            if(cond){
              if(smooth == TRUE){
                sm.df <- smoothEachPosition(features = object@features,
                                            posdf = ribo.df,
                                            slide_window = slide_window)
              }else{
                sm.df <- ribo.df
              }
            }else{
              sm.df <- ribo.df
            }

            if(smooth == TRUE){
              if(serp_exp == "total"){
                object@total_occupancy <- sm.df
              }else{
                object@ip_occupancy <- sm.df
              }
            }else{
              sm.df$smooth <- sm.df$rpm

              if(serp_exp == "total"){
                object@total_occupancy <- sm.df
              }else{
                object@ip_occupancy <- sm.df
              }
            }

            object@gene_name <- gene_name

            return(object)
          }
)



