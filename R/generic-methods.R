# ==============================================================================
# function to get detailed information for reads
# ==============================================================================

#' Generate Summary Information for Ribosome Profiling Data
#'
#' @description Generic function to generate summary information from ribosome profiling
#' or RNA-seq data.
#'
#' @param object An object containing experimental data
#' @param ... Additional arguments passed to methods
#'
#' @return An updated object containing summary information
#'
#' @export
#' @rdname generate_summary
setGeneric("generate_summary",function(object,...) standardGeneric("generate_summary"))



#' Generate summary statistics for ribosome profiling or RNA-seq data
#'
#' @description
#' This method processes BAM files to generate position-wise read counts and mapping statistics
#' for either ribosome profiling or RNA-seq data. It handles both transcriptome and genome-aligned
#' data, performing appropriate coordinate conversions when necessary.
#'
#' @param object A ribotrans object containing BAM file information and transcript features
#' @param exp_type Character, specify the experiment type to analyze:
#'   - "ribo": ribosome profiling data
#'   - "rna": RNA-seq data
#' @param nThreads Integer, number of threads to use for BAM file processing
#'
#' @details
#' The function performs the following steps:
#' 1. Filters library information based on experiment type
#' 2. Processes each BAM file to extract read positions and counts
#' 3. For genome-aligned data:
#'    - Converts genomic coordinates to transcript coordinates
#'    - Handles 5' or 3' end assignment based on strand information
#' 4. For transcriptome-aligned data:
#'    - Directly uses transcript coordinates
#' 5. Merges position information with transcript annotations
#'
#' The resulting summary includes:
#' - rname: transcript identifier
#' - pos: position along transcript
#' - qwidth: read length
#' - count: number of reads at each position
#' - sample: sample identifier
#' - mstart/mstop: CDS start/stop positions
#' - translen: transcript length
#'
#' @return Returns the input ribotrans object with updated summary_info slot containing
#' the processed mapping statistics
#'
#' @importFrom Rsamtools scanBam ScanBamParam scanBamFlag
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges findOverlaps
#' @importFrom dplyr mutate rename group_by summarise n case_when
#' @importFrom fastplyr f_group_by f_summarise f_left_join
#' @importFrom purrr map_df
#'
#' @note
#' - For large BAM files, increasing nThreads can significantly improve processing speed
#' - Memory usage increases with BAM file size and number of threads
#' - Genome-aligned data requires more processing time due to coordinate conversion
#'
#' @seealso
#' \code{\link{construct_ribotrans}} for creating ribotrans objects
#' \code{\link{get_occupancy}} for calculating ribosome occupancy
#'
#' @export
#' @rdname generate_summary
setMethod("generate_summary",
          signature(object = "ribotrans"),
          function(object, exp_type = c("ribo", "rna"),
                   nThreads = 1){
            exp_type <- match.arg(exp_type, c("ribo", "rna"))

            bams <- subset(object@library, type == exp_type)
            bfn <- bams$bam
            gp <- bams$sample

            features <- object@features

            # x = 1
            purrr::map_df(seq_along(bfn),function(x){
              # get position
              bam_data <- Rsamtools::scanBam(file = bfn[x],
                                             nThreads = nThreads,
                                             param = Rsamtools::ScanBamParam(what = c("rname", "pos", "strand", "qwidth"),
                                                                             flag = Rsamtools::scanBamFlag(isUnmappedQuery = FALSE))
              )

              tinfo <- data.frame(bam_data[[1]])

              # check mapping_type
              if(object@mapping_type == "genome"){
                if(object@assignment_mode == "end5"){
                  # check strand to assign 5'end position
                  tinfo <- tinfo %>%
                    dplyr::mutate(pos = ifelse(strand == "-", pos + qwidth - 1, pos))
                }else{
                  # check strand to assign 3'end position
                  tinfo <- tinfo %>%
                    dplyr::mutate(pos = ifelse(strand == "+", pos + qwidth - 1, pos))
                }

                tinfo.gr <- tinfo %>%
                  fastplyr::f_group_by(rname,pos,qwidth) %>%
                  fastplyr::f_summarise(count = dplyr::n(),.groups = "drop") %>%
                  dplyr::rename(seqnames = rname, start = pos) %>%
                  dplyr::mutate(end = start,.after = "start")

                # to GRanges
                tinfo.gr <- GenomicRanges::GRanges(tinfo.gr)

                # transcript annotation regions
                trans_rg <- GenomicRanges::GRanges(object@genome_trans_features)

                # genomic position to transcriptome position
                ov <- IRanges::findOverlaps(query = tinfo.gr,subject = trans_rg)

                # get overlap data
                if (requireNamespace("S4Vectors", quietly = TRUE)) {
                  lo <- cbind(as.data.frame(tinfo.gr[S4Vectors::queryHits(ov)]),
                              as.data.frame(trans_rg[S4Vectors::subjectHits(ov)]))
                } else {
                  warning("Package 'S4Vectors' is needed for this function to work.")
                }


                # make unique names
                names(lo) <- make.names(names(lo),unique = T)

                # calculate transcript position
                tinfo <- lo %>%
                  dplyr::mutate(pos = dplyr::case_when(strand.1 == "+" ~ tx_len - abs(end.1 - start),
                                                       strand.1 == "-" ~ tx_len - abs(start - start.1)),
                                rname = paste(transcript_id,gene_name,sep = "|")) %>%
                  fastplyr::f_select(rname,pos,qwidth,count)

              }else{
                tinfo <- tinfo %>%
                  fastplyr::f_group_by(rname,pos,qwidth) %>%
                  fastplyr::f_summarise(count = dplyr::n())
              }

              tinfo$sample <- gp[x]

              return(tinfo)
            }) -> all.tinfo

            # merge with transcript annotation
            tinfo.anno <- all.tinfo %>%
              fastplyr::f_left_join(y = features[,c("idnew","mstart","mstop","translen")],
                                    by = c("rname" = "idnew"))


            object@summary_info <- tinfo.anno
            return(object)
          }
)





# ==============================================================================
# function to calculate coverage
# ==============================================================================

#' @title Get RNA-seq Coverage from a `ribotrans` Object
#' @description This generic function extracts RNA-seq coverage from a `ribotrans` object.
#' @param object A `ribotrans` object containing RNA-seq data.
#' @param ... Additional arguments passed to methods.
#' @return A modified `ribotrans` object with RNA-seq coverage information.
#' @seealso \code{\link{get_coverage,ribotrans-method}} for the specific implementation.
#' @export
#' @rdname get_coverage
setGeneric("get_coverage",function(object,...) standardGeneric("get_coverage"))



#' @title Extract RNA Coverage for a `ribotrans` Object
#' @description This method extracts RNA-seq coverage information from a `ribotrans` object,
#'     optionally applying a smoothing transformation.
#'
#' @param object A `ribotrans` object containing RNA-seq BAM files and library metadata.
#' @param gene_name Character. Name of the gene from which to extract coverage.
#'    If `NULL` (default), retrieves coverage for all genes.
#' @param smooth Character. Should the data be smoothed? Options: `FALSE` or `TRUE`. Default is `FALSE`.
#' If coordinate_to_trans is FALSE and mapping_type is genome, data will not be smoothed.
#' @param coordinate_to_trans Logical. Whether to convert genomic coordinates to transcript
#'        coordinates. Default is FALSE.
#' @param slide_window Integer. The window size for smoothing (only used if `smooth = "TRUE"`). Default: `30`.
#'
#' @details
#' RNA-seq coverage describes read distribution along transcripts, providing insights into gene expression levels.
#' This function processes `ribotrans` BAM files and computes RNA read coverage for each transcript.
#'
#' **Analysis Steps:**
#' 1. Extracts RNA-seq BAM files from `object@bam_file`, filtering by `"rna"` type.
#' 2. Retrieves mapped read counts from `object@library` to normalize coverage.
#' 3. Calls:
#'    - `getCoverageGenome()` for **genome-mapped BAM files**.
#'    - `getCoverage()` for **transcript-mapped BAM files**.
#' 4. Merges all samples into a single `data.frame`.
#' 5. **If `smooth = "TRUE"`**, applies a **sliding window smoothing** (`smoothEachPosition()`).
#'
#' The processed data is stored in `object@RNA_coverage`, and the smoothing status is recorded in `object@RNA.smoothed`.
#'
#' @return A modified `ribotrans` object that now contains:
#' - **`object@RNA_coverage`**: A `data.frame` with RNA-seq coverage values.
#' - **`object@RNA.smoothed`**: `"TRUE"` or `"FALSE"`, indicating whether coverage data was smoothed.
#'
#' @seealso
#' - \code{\link{getCoverageGenome}}: Extract RNA coverage from genome-mapped data.
#' - \code{\link{getCoverage}}: Extract RNA coverage from transcript-aligned data.
#' - \code{\link{ribotrans}}: The S4 object storing this data.
#'
#' @examples
#' \dontrun{
#' # Assume ribo_obj is a ribotrans object with RNA BAM files
#' result <- get_coverage(ribo_obj, gene_name = "TP53", smooth = "TRUE", slide_window = 50)
#' print(result@RNA_coverage)
#' }
#'
#' @importFrom purrr map_df
#' @importFrom dplyr select filter mutate
#' @export
#' @rdname get_coverage
setMethod("get_coverage",
          signature(object = "ribotrans"),
          function(object, gene_name = NULL,
                   smooth = FALSE,
                   coordinate_to_trans = FALSE,
                   slide_window = 30){
            # smooth <- match.arg(smooth, c(FALSE, TRUE))

            bf <- subset(object@bam_file, type == "rna")
            rnabams <- bf$bam
            gp <- bf$sample

            lib <- subset(object@library, type == "rna")
            lib <- lib[match(rnabams, lib$bam),]
            totalreads <- lib$mappped_reads

            # loop for each sample
            purrr::map_df(seq_along(rnabams),function(x){
              # check mapping type
              if(object@mapping_type == "genome"){
                tmp <- getCoverageGenome(bam_file = rnabams[x],
                                         gene_name = gene_name,
                                         features = object@genome_trans_features,
                                         total_reads = totalreads[x],
                                         coordinate_to_trans = coordinate_to_trans)
              }else{
                tmp <- getCoverage(bam_file = rnabams[x],
                                   gene_name = gene_name,
                                   features = object@features,
                                   total_reads = totalreads[x])
              }

              # add sample name
              tmp$sample <- gp[x]

              return(tmp)
            }) -> rna.df

            # smooth for each position
            if(coordinate_to_trans == TRUE){
              sm.df <- smoothEachPosition(features = object@features,
                                          posdf = rna.df,
                                          slide_window = slide_window)
            }else{
              sm.df <- rna.df
            }

            if(smooth == TRUE){
              object@RNA_coverage <- sm.df
              object@RNA.smoothed <- "TRUE"
            }else{
              sm.df$smooth <- sm.df$rpm
              object@RNA_coverage <- sm.df
              object@RNA.smoothed <- "FALSE"
            }

            object@gene_name <- gene_name

            return(object)
          }
)



# ==============================================================================
# function to calculate ribosome occupancy
# ==============================================================================

#' @title Get Ribosome Occupancy from `ribotrans` Object
#' @description This generic function retrieves ribosome occupancy from a `ribotrans` object.
#' @param object A `ribotrans` object containing ribosome profiling data.
#' @param ... Additional arguments passed to specific methods.
#' @return A modified `ribotrans` object with ribosome occupancy information.
#' @seealso \code{\link{get_occupancy,ribotrans-method}} for the specific method.
#' @export
#' @rdname get_occupancy
setGeneric("get_occupancy",function(object,...) standardGeneric("get_occupancy"))


#' @title Extract Ribosome Occupancy for a `ribotrans` Object
#' @description This method extracts ribosome occupancy information from a `ribotrans` object,
#'     optionally applying a smoothing transformation.
#'
#' @param object A `ribotrans` object containing ribosome profiling BAM files and library metadata.
#' @param gene_name Character. Name of the gene from which to extract coverage.
#'    If `NULL` (default), retrieves coverage for all genes.
#' @param smooth Character. Should the data be smoothed? Options: `FALSE` or `TRUE`. Default is `FALSE`.
#' If coordinate_to_trans is FALSE and mapping_type is genome, data will not be smoothed.
#' @param coordinate_to_trans Logical. Whether to convert genomic coordinates to transcript
#'        coordinates. Default is FALSE.
#' @param do_reads_offset Logical; if set to \code{TRUE}, adjusts read positions based on
#' precomputed offset information stored in the \code{reads_offset_info} slot (default: \code{FALSE}).
#' @param slide_window Integer. The number of nucleotides for smoothing (used only if `smooth = "TRUE"`). Default is `30`.
#'
#' @details
#' Ribosome occupancy describes the density of ribosomes along a transcript.
#' This function processes `ribotrans` BAM files to compute read coverage, supporting both **genome-aligned**
#' and **transcriptome-aligned** data.
#'
#' **Analysis Steps:**
#' 1. Extracts mapped ribosome profiling BAM files from `object@bam_file` (filtered by `"ribo"` type).
#' 2. Retrieves mapped read counts from `object@library` to normalize the occupancy.
#' 3. Calls:
#'    - `getOccupancyGenome()` for **genome-based BAM files**.
#'    - `getOccupancy()` for **transcript-based BAM files**.
#' 4. Merges all samples into a single `data.frame`.
#' 5. **If `smooth = "TRUE"`**, applies a **sliding window smoothing** (`smoothEachPosition()`).
#'
#'
#' The processed data is stored in `object@ribo_occupancy`, and the smoothing status is recorded in `object@ribo.smoothed`.
#'
#' @return A modified `ribotrans` object that now contains:
#' - **`object@ribo_occupancy`**: A `data.frame` with ribosome occupancy values.
#' - **`object@ribo.smoothed`**: `"TRUE"` or `"FALSE"`, indicating whether occupancy data was smoothed.
#'
#' @seealso
#' - \code{\link{getOccupancyGenome}}: Extract occupancy from genome-mapped data.
#' - \code{\link{getOccupancy}}: Extract occupancy from transcript-aligned data.
#' - \code{\link{ribotrans}}: The S4 object storing this data.
#'
#' @examples
#' \dontrun{
#' # Assume ribo_obj is a ribotrans object with ribosome BAM files
#' result <- get_occupancy(ribo_obj, gene_name = "TP53", smooth = "TRUE", slide_window = 50)
#' print(result@ribo_occupancy)
#' }
#'
#' @importFrom purrr map_df
#' @importFrom dplyr select filter mutate
#' @export
#' @rdname get_occupancy
setMethod("get_occupancy",
          signature(object = "ribotrans"),
          function(object, gene_name = NULL,
                   smooth = FALSE,
                   coordinate_to_trans = FALSE,
                   do_reads_offset = FALSE,
                   slide_window = 30){
            # smooth <- match.arg(smooth, c(FALSE, TRUE))
            assignment_mode <- object@assignment_mode

            bf <- subset(object@bam_file, type == "ribo")
            ribobams <- bf$bam
            gp <- bf$sample

            lib <- subset(object@library, type == "ribo")
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
            if(coordinate_to_trans == TRUE){
              sm.df <- smoothEachPosition(features = object@features,
                                          posdf = ribo.df,
                                          slide_window = slide_window)
            }else{
              sm.df <- ribo.df
            }


            if(smooth == TRUE){
              object@ribo_occupancy <- sm.df
              object@ribo.smoothed <- "TRUE"
            }else{
              sm.df$smooth <- sm.df$rpm
              object@ribo_occupancy <- sm.df
              object@ribo.smoothed <- "FALSE"
            }

            object@gene_name <- gene_name

            return(object)
          }
)



# ==============================================================================
# scale to mRNA abundance
# ==============================================================================

#' @title Compute Scaled Ribosome Occupancy
#' @description This generic function calculates scaled ribosome occupancy by normalizing ribosome occupancy
#'     with RNA-seq coverage.
#' @param object A `ribotrans` object containing RNA-seq coverage and ribosome occupancy.
#' @param ... Additional arguments passed to specific methods.
#' @return A modified `ribotrans` object with scaled occupancy stored in `object@scaled_occupancy`.
#' @seealso \code{\link{get_scaled_occupancy,ribotrans-method}} for the method implementation.
#' @export
#' @rdname get_scaled_occupancy
setGeneric("get_scaled_occupancy",function(object,...) standardGeneric("get_scaled_occupancy"))


#' @title Compute Scaled Ribosome Occupancy for a `ribotrans` Object
#' @description This method calculates the scaled ribosome occupancy by dividing ribosome density
#'     (`ribo_occupancy`) by RNA-seq coverage (`RNA_coverage`), optionally applying a smoothing function.
#'
#' @param object A `ribotrans` object containing both `ribo_occupancy` and `RNA_coverage` data.
#' @param smooth Character. Should the data be smoothed? Options: `TRUE` or `FALSE`. Default is `TRUE`.
#' @param slide_window Integer. The window size for smoothing (used only if `smooth = "TRUE"`). Default: `30`.
#'
#' @details
#' **Functionality:**
#' This method normalizes ribosome occupancy by RNA-seq coverage per position, producing a **scaled occupancy profile**.
#'
#' **Workflow:**
#' 1. Ensure both `RNA_coverage` (`object@RNA_coverage`) and `ribo_occupancy` (`object@ribo_occupancy`) exist.
#' 2. Compute: `scaled occupancy = (ribo occupancy) / (RNA coverage)`.
#' 3. Replace `NA` or `Inf` values with `0` for stability.
#' 4. **Optional:** Apply a rolling mean (`zoo::rollmean()`) for smoothing.
#'
#' **Error Handling:**
#' - If `RNA_coverage` or `ribo_occupancy` is missing, the function stops with an error.
#'
#' @return A modified `ribotrans` object containing `object@scaled_occupancy`, a `data.frame` with the scaled enrichment values.
#' @seealso
#' \code{\link{get_occupancy}}: Extract raw ribosome occupancy.
#' \code{\link{get_coverage}}: Extract RNA coverage.
#' \code{\link{ribotrans}}: The S4 class storing input data.
#'
#' @examples
#' \dontrun{
#' # Assume `ribo_obj` is a ribotrans object with RNA and ribosome data
#' result <- get_scaled_occupancy(ribo_obj, smooth = "TRUE", slide_window = 50)
#' print(result@scaled_occupancy)
#' }
#'
#' @importFrom dplyr select left_join mutate filter
#' @importFrom purrr map_df
#' @export
#' @rdname get_scaled_occupancy
setMethod("get_scaled_occupancy",
          signature(object = "ribotrans"),
          function(object,
                   smooth = TRUE,
                   slide_window = 30){
            # smooth <- match.arg(smooth, c(TRUE, FALSE))

            # ribo rpm/rna rpm for each position

            # check rna coverage
            if(is.null(object@RNA_coverage)){
              stop("Please supply RNA coverage data!")
            }

            if(is.null(object@ribo_occupancy)){
              stop("Please supply Ribo occupancy data!")
            }

            # run
            mb <- object@ribo_occupancy %>%
              dplyr::select(sample,rname,pos,rpm) %>%
              dplyr::left_join(y = object@RNA_coverage[,c("sample","rname","pos","rpm")],
                               by = c("sample","rname","pos")) %>%
              dplyr::mutate(enrich = (rpm.x/rpm.y)) %>%
              dplyr::mutate(enrich = ifelse(is.na(enrich) | is.infinite(enrich),0,enrich))

            # loop for each sample and transcript to smooth
            tanno <- subset(object@features, gene == object@gene_name)

            # check whether smooth data
            if(smooth == TRUE){
              purrr::map_df(1:nrow(tanno),function(x){
                tmp <- tanno[x,]

                tmp2 <- mb %>%
                  dplyr::filter(rname == tmp$idnew)

                sp <- unique(tmp2$sample)

                # loop for each sample
                purrr::map_df(seq_along(sp),function(x){
                  tmp3 <- tmp2 %>%
                    dplyr::filter(sample == sp[x])

                  # smooth data
                  if (requireNamespace("zoo", quietly = TRUE)) {
                    tmp3$smooth <- zoo::rollmean(tmp3$enrich, k = slide_window, fill = 0)
                  } else {
                    warning("Package 'zoo' is needed for this function to work.")
                  }


                  return(tmp3)
                }) -> smooth.df

                return(smooth.df)
              }) -> mb.smoothed
            }else{
              mb$smooth <- mb$enrich
              mb.smoothed <- mb
            }


            object@scaled_occupancy <- mb.smoothed
            return(object)
          }
)




