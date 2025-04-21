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
#' @param slide_window Numeric. Size of the sliding window for smoothing. Default: 30.
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
                   slide_window = 30){
            serp_exp <- match.arg(serp_exp,choices = c("total","ip"))

            assignment_mode <- object@assignment_mode

            bf <- subset(object@bam_file, type == serp_exp)
            ribobams <- bf$bam
            gp <- bf$sample

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
              tmp$sample <- gp[x]

              # ================================================================
              # whether do read offset
              if(do_reads_offset == TRUE){
                # check offset data
                if(nrow(object@reads_offset_info) == 0){
                  stop("No reads offset information in object!")
                }

                offset <- object@reads_offset_info

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

              }else{
                adjusted <- tmp
              }

              adjusted <- adjusted %>%
                dplyr::group_by(sample,rname,pos) %>%
                dplyr::summarise(count = sum(count),rpm = sum(rpm), .groups = "drop")

              return(adjusted)
            }) -> ribo.df

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



# ==============================================================================
# do enrichment analysis
# ==============================================================================

#' Calculate Enrichment from Ribosome Profiling Data
#'
#' @description
#' Generic function for calculating enrichment ratios between IP and total samples in selective
#' ribosome profiling experiments.
#'
#' @param object An S4 object containing ribosome profiling data.
#' @param ... Additional arguments passed to specific method implementations.
#'
#' @return An updated object with calculated enrichment data.
#'
#' @details
#' This is a generic function that computes enrichment ratios between immunoprecipitated (IP)
#' and total ribosome profiling samples. The specific implementation depends on the class of the
#' input object.
#'
#' @examples
#' \dontrun{
#' # Create a serp object and process data
#' serp_obj <- construct_serp(...)
#' serp_obj <- get_occupancy(serp_obj, serp_exp = "total")
#' serp_obj <- get_occupancy(serp_obj, serp_exp = "ip")
#'
#' # Calculate enrichment
#' serp_obj <- get_enrichment(serp_obj)
#' }
#'
#' @export
setGeneric("get_enrichment",function(object,...) standardGeneric("get_enrichment"))



#' @describeIn get_enrichment Calculate enrichment ratios for SERP objects
#'
#' @param smooth Logical. Whether to apply smoothing to the enrichment profile. Default: TRUE.
#' @param window_size Numeric. Size of the sliding window for smoothing. Default: 105.
#' @param ... Additional arguments (currently unused).
#'
#' @details
#' For objects of class 'serp', this method calculates the enrichment ratio between
#' IP and total ribosome profiling samples. The method requires both 'total_occupancy'
#' and 'ip_occupancy' data to be present in the object.
#'
#' The enrichment calculation process includes:
#' 1. Calculating the ratio of IP RPM to total RPM for each position
#' 2. Filtering for coding sequence (CDS) regions
#' 3. Optionally applying smoothing using a sliding window approach
#' 4. When smoothing is applied, confidence intervals are calculated using the Agresti-Coull method
#'    and normalized by the ratio of total mapped reads between IP and total samples
#' 5. When smoothing is not applied, positions are grouped by codon (every 3 nucleotides)
#'    and average enrichment is calculated
#'
#' @return An updated 'serp' object with enrichment data stored in the 'enriched_ratio' slot.
#'
#' @examples
#' \dontrun{
#' # Process a serp object
#' serp_obj <- get_occupancy(serp_obj, serp_exp = "total")
#' serp_obj <- get_occupancy(serp_obj, serp_exp = "ip")
#'
#' # Calculate enrichment with default smoothing
#' serp_obj <- get_enrichment(serp_obj)
#'
#' # Calculate enrichment without smoothing
#' serp_obj <- get_enrichment(serp_obj, smooth = FALSE)
#'
#' # Calculate enrichment with custom window size
#' serp_obj <- get_enrichment(serp_obj, window_size = 51)
#' }
#'
#' @importFrom dplyr select left_join mutate filter group_by summarise
#' @importFrom purrr map_df
setMethod("get_enrichment",
          signature(object = "serp"),
          function(object,
                   smooth = TRUE,
                   window_size = 105){

            # ribo rpm/rna rpm for each position

            # check rna coverage
            if(is.null(object@total_occupancy)){
              stop("Please supply total_occupancy data!")
            }

            if(is.null(object@ip_occupancy)){
              stop("Please supply ip_occupancy data!")
            }

            # run
            mb <- object@total_occupancy %>%
              dplyr::select(sample,rname,pos,count,rpm) %>%
              dplyr::left_join(y = object@ip_occupancy[,c("sample","rname","pos","count","rpm")],
                               by = c("sample","rname","pos")) %>%
              dplyr::mutate(enrich = (rpm.y/rpm.x)) %>%
              dplyr::mutate(enrich = ifelse(is.na(enrich) | is.infinite(enrich),0,enrich))

            mb[is.na(mb)] <- 0

            # loop for each sample and transcript to smooth
            tanno <- subset(object@features, gene == object@gene_name)

            # check whether smooth data
            # x = 1
            purrr::map_df(1:nrow(tanno),function(x){
              tmp <- tanno[x,]

              tmp2 <- mb %>%
                dplyr::filter(rname == tmp$idnew)

              sp <- unique(tmp2$sample)

              # loop for each sample
              # x = 1
              purrr::map_df(seq_along(sp),function(x){
                tmp3 <- tmp2 %>%
                  dplyr::filter(sample == sp[x])

                # each pos
                pos.df <- data.frame(rname = tmp3$rname[1],pos = 1:tmp$exonlen)

                # merge with value
                pos.df <- pos.df %>% dplyr::left_join(y = tmp3,by = c("rname", "pos"))
                pos.df$sample <- tmp3$sample[1]
                pos.df[is.na(pos.df)] <- 0

                # retain CDS region
                pos.df <- pos.df %>%
                  dplyr::filter(pos >= tmp$mstart & pos <= tmp$mstop) %>%
                  dplyr::mutate(pos = pos - tmp$mstart + 1)

                # ==============================================================
                if(smooth == TRUE){
                  # smooth data for each position
                  ip_total <- subset(object@library, sample == sp[x] & type == "ip")$mappped_reads
                  total_total <- subset(object@library, sample == sp[x] & type == "total")$mappped_reads

                  probnorm <- ip_total / total_total

                  len_gene <- tmp$cds
                  ip <- pos.df$count.y
                  total <- pos.df$count.x

                  win_starts <- 1:(len_gene - window_size + 1)
                  win_ends <- window_size:len_gene

                  cumsum_ip <- c(0, cumsum(ip))
                  win_ip <- cumsum_ip[win_ends+1] - cumsum_ip[win_starts]

                  cumsum_total <- c(0, cumsum(total))
                  win_total <- cumsum_total[win_ends+1] - cumsum_total[win_starts]

                  if (requireNamespace("binom", quietly = TRUE)) {
                    cidf <- binom::binom.agresti.coull(win_ip, win_ip + win_total, conf.level = 0.95)
                  } else {
                    warning("Package 'binom' is needed for this function to work.")
                  }

                  lower <- prob2odds(pmax(0, cidf$lower))
                  lower <- ifelse(is.nan(lower), 0, lower)

                  upper <- prob2odds(pmin(1, cidf$upper))
                  upper <- ifelse(is.nan(upper), Inf, upper)

                  mean <- prob2odds(cidf$mean)
                  mean <- ifelse(is.nan(mean), 0, mean)

                  # output
                  df <- data.frame(sample = sp[x],
                                   rname = pos.df$rname[x],
                                   pos = win_starts,
                                   win_ip = win_ip,
                                   win_total = win_total,
                                   enrich = mean/probnorm,
                                   lower = lower/probnorm,
                                   upper = upper/probnorm)
                }else{
                  df <- pos.df %>%
                    dplyr::select(sample, rname, pos, enrich) %>%
                    dplyr::mutate(pos = (pos - 1) %/% 3 + 1) %>%
                    dplyr::group_by(sample,rname,pos) %>%
                    dplyr::summarise(enrich = mean(enrich))
                }

                return(df)
              }) -> smooth.df

              return(smooth.df)
            }) -> mb.smoothed

            mb.smoothed$window_size <- window_size

            object@enriched_ratio <- mb.smoothed
            return(object)
          }
)
