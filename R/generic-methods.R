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



#' generate_summary Generate summary for ribotrans objects
#'
#' @description This method extracts read mapping information from BAM files and links it
#' to transcript features in a `ribotrans` object.
#'
#' @param object A `ribotrans` object containing library and feature information.
#' @param exp_type Character string indicating experiment type. One of `"ribo"` (ribosome profiling)
#' or `"rna"` (RNA-seq).
#' @param nThreads Integer specifying the number of threads to use for BAM file processing.
#'
#' @return The input `ribotrans` object with an updated `summary_info` slot containing read count summaries.
#'
#' @details This method processes BAM files corresponding to the specified experiment type (`ribo` or `rna`) and extracts:
#' \itemize{
#'   \item Read positions (`pos`)
#'   \item Read lengths (`qwidth`)
#'   \item Read counts per transcript position
#' }
#' The extracted mapping data is then merged with transcript annotation data from the `features` slot.
#'
#' @importFrom Rsamtools scanBam ScanBamParam scanBamFlag
#' @importFrom purrr map_df
#' @importFrom fastplyr f_group_by f_summarise f_left_join
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
                                             param = Rsamtools::ScanBamParam(what = c("rname", "pos", "qwidth"),
                                                                             flag = Rsamtools::scanBamFlag(isUnmappedQuery = FALSE))
              )

              tinfo <- data.frame(bam_data[[1]]) %>%
                fastplyr::f_group_by(rname,pos,qwidth) %>%
                fastplyr::f_summarise(count = n())

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
#' @param smooth Character. Should the data be smoothed? Options: `"FALSE"` or `"TRUE"`. Default is `"FALSE"`.
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
                   smooth = c("FALSE", "TRUE"),
                   slide_window = 30){
            smooth <- match.arg(smooth, c("FALSE", "TRUE"))

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
                                         total_reads = totalreads[x])
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
            if(smooth == TRUE){
              sm.df <- smoothEachPosition(features = object@features,
                                          posdf = rna.df,
                                          slide_window = slide_window)

              object@RNA_coverage <- sm.df
              object@RNA.smoothed <- "TRUE"
            }else{
              rna.df$smooth <- rna.df$rpm
              object@RNA_coverage <- rna.df
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
#' @param smooth Character. Should the data be smoothed? Options: `"FALSE"` or `"TRUE"`. Default is `"FALSE"`.
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
                   smooth = c("FALSE", "TRUE"),
                   slide_window = 30){
            smooth <- match.arg(smooth, c("FALSE", "TRUE"))

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
                                          total_reads = totalreads[x])
              }else{
                tmp <- getOccupancy(bam_file = ribobams[x],
                                    gene_name = gene_name,
                                    features = object@features,
                                    total_reads = totalreads[x])
              }

              # add sample name
              tmp$sample <- gp[x]

              return(tmp)
            }) -> ribo.df

            # smooth for each position
            if(smooth == TRUE){
              sm.df <- smoothEachPosition(features = object@features,
                                          posdf = ribo.df,
                                          slide_window = slide_window)

              object@ribo_occupancy <- sm.df
              object@ribo.smoothed <- "TRUE"
            }else{
              ribo.df$smooth <- ribo.df$rpm
              object@ribo_occupancy <- ribo.df
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
#' @param smooth Character. Should the data be smoothed? Options: `"TRUE"` or `"FALSE"`. Default is `"TRUE"`.
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
#' @importFrom zoo rollmean
#' @export
#' @rdname get_scaled_occupancy
setMethod("get_scaled_occupancy",
          signature(object = "ribotrans"),
          function(object,
                   smooth = c("TRUE", "FALSE"),
                   slide_window = 30){
            smooth <- match.arg(smooth, c("TRUE", "FALSE"))

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
              dplyr::select(sample,rname,pos,smooth) %>%
              dplyr::left_join(y = object@RNA_coverage[,c("sample","rname","pos","smooth")],
                               by = c("sample","rname","pos")) %>%
              dplyr::mutate(enrich = (smooth.x/smooth.y)) %>%
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
                  tmp3$smooth <- zoo::rollmean(tmp3$enrich, k = slide_window, fill = 0)

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




