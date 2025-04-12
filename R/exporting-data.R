# ==============================================================================
# function to export bigwig data for ribo occupancy from bams on genome
# ==============================================================================

#' @title Export Genome-Wide Ribosome Occupancy
#'
#' @description
#' `export_genome_occupancy()` is a generic function designed to export genome-wide
#' ribosome footprint occupancy data in **BigWig format (BW)**. It supports **E/P/A-site correction**
#' and normalization for ribosome profiling data.
#'
#' @param object An object of class **`ribotrans`** containing Ribo-Seq data.
#' @param ... Additional arguments passed to the specific method.
#'
#' @details
#' This function is a **generic function**, and its specific implementation depends
#' on the **class of `object`**. See **methods listed under `showMethods("export_genome_occupancy")`**.
#'
#' The `export_genome_occupancy()` method is mainly used for exporting
#' **ribosome occupancy signals** from BAM files to **BigWig format (BW)**
#' for visualization in tools like **IGV or UCSC Genome Browser**.
#'
#' @return
#' The specific method implementation for **`ribotrans` objects** exports **BigWig files** (BW)
#' to the specified directory and **returns a confirmation message**.
#'
#' @seealso
#' - `export_genome_occupancy,ribotrans-method` for the specific `ribotrans` method.
#' - `getOccupancyGenome()` for extracting genome-wide ribosome occupancy.
#'
#' @export
setGeneric("export_genome_occupancy",function(object,...) standardGeneric("export_genome_occupancy"))




#' @title Export Ribosome Footprint Data as BigWig
#'
#' @description
#' This method extracts **ribosome footprint occupancy data from BAM files**
#' and exports it as **BigWig (.bw) format**, supporting **P-site offset correction**
#' and **RPM normalization**.
#'
#' @param object A **`ribotrans`** object containing Ribo-Seq experimental data.
#' @param genome_file A **character string** specifying the path to a **FASTA file**
#' containing reference genome sequences.
#' @param output_path A **character string** specifying the output directory for BigWig files
#' (default: `"./"`).
#' @param do_reads_offset Logical; if `TRUE`, applies **E/P/A-site offset correction**
#' using metadata in `object@reads_offset_info` (default: `FALSE`).
#' @param shift Integer defining the **offset shift** applied to read positions. **Default**: `0`.
#'
#' @details
#' This method reads Ribo-Seq data from BAM files, extracts **read occupancy**,
#' applies optional **P-site offset correction**, and **normalizes data** in **RPM (Reads Per Million)**.
#' Finally, it exports the processed data as **BigWig (.bw) format** for genome-wide visualization.
#'
#' **Steps:**
#' 1. Extract **mapped reads** from BAM files.
#' 2. Assign **E/P/A-site positions** based on strand and fragment length.
#' 3. Apply **normalization** and **offset correction** (if enabled).
#' 4. Save data in **BigWig (bw) format**.
#'
#' ⚠️ **Note:** If `object@mapping_type != "transcriptome"`, an error is returned.
#'
#' @return
#' - **BigWig files (.bw)** are saved in the `output_path` directory.
#' - A **confirmation message** is printed upon completion.
#'
#' @examples
#' \dontrun{
#' ribo_obj <- load_ribotrans("ribo_data.rds")
#'
#' # Export BigWig files (without offset correction)
#' export_genome_occupancy(ribo_obj, genome_file = "reference.fa",
#'                         output_path = "bw_output/", do_reads_offset = FALSE)
#'
#' # Apply E/P/A-site offset correction before exporting
#' export_genome_occupancy(ribo_obj, genome_file = "reference.fa",
#'                         output_path = "bw_output/", do_reads_offset = TRUE)
#' }
#'
#' @seealso
#' - `getOccupancyGenome()` for extracting genome-wide ribosome occupancy.
#' - `rtracklayer::export.bw()` for exporting BigWig files.
#'
#' @note Requires `rtracklayer`, `purrr`, and `GenomicRanges` packages.
#'
#' @export
setMethod("export_genome_occupancy",
          signature(object = "ribotrans"),
          function(object,
                   genome_file = NULL,
                   output_path = "./",
                   do_reads_offset = FALSE,
                   shift = 0){
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
                                          total_reads = totalreads[x],
                                          assignment_mode = assignment_mode,
                                          coordinate_to_trans = FALSE,
                                          return_all_position = TRUE)
              }else{
                warnings("mapping_type for transcriptome is not supoorted!")
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

                adjusted <- tmp %>%
                  dplyr::left_join(y = offset[,c("sample","rel_pos","qwidth")],by = c("sample","qwidth")) %>%
                  stats::na.omit() %>%
                  dplyr::mutate(pos = dplyr::if_else(strand == "+",
                                                     pos - rel_pos + shift,pos + rel_pos - shift)) %>%
                  dplyr::select(-strand,-qwidth,-rel_pos)

              }else{
                adjusted <- tmp
              }

              # to GRanges
              adjusted <- adjusted %>%
                fastplyr::f_group_by(sample,rname,pos) %>%
                fastplyr::f_summarise(count = sum(count),rpm = sum(rpm), .groups = "drop") %>%
                fastplyr::f_select(rname,pos,rpm) %>%
                dplyr::rename(seqnames = rname, start = pos,score = rpm) %>%
                dplyr::mutate(end = start,.after = "start") %>%
                GenomicRanges::GRanges()

              # ================================================================
              # extract chromosome length
              fa <- Biostrings::readDNAStringSet(filepath = genome_file)
              names(fa) <- sapply(strsplit(names(fa),split = " "),"[",1)


              # filter seqnames in target
              seqlengths.info <- data.frame(chr = names(GenomeInfoDb::seqlengths(fa)),
                                            len = GenomeInfoDb::seqlengths(fa)) %>%
                dplyr::filter(chr %in% names(GenomeInfoDb::seqlengths(adjusted)))

              # order
              seqlengths.info <- seqlengths.info[names(GenomeInfoDb::seqlengths(adjusted)),]


              # output
              if (requireNamespace("rtracklayer", quietly = TRUE)) {
                # assign seqlengths information for adjusted
                GenomeInfoDb::seqlengths(adjusted) <- seqlengths.info$len

                rtracklayer::export.bw(object = adjusted,con = paste(output_path,gp[x],"_occuapcny.bw",sep = ""))
              } else {
                warning("Package 'rtracklayer' is needed for this function to work.")
              }


              return(message(paste(gp[x],"has been finished.")))
            })


          }
)
