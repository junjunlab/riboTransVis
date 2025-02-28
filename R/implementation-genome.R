
#' Extract Transcript Information from a GTF File
#'
#' This function loads a GTF file and extracts exon information for each transcript.
#' It computes the total transcript length and determines the exon positions
#' relative to the transcript.
#'
#' @param gtf_file Character. Path to the GTF file.
#'
#' @return A `data.frame` containing transcript position information with the columns:
#' \itemize{
#'   \item `seqnames`: Chromosome name.
#'   \item `start`, `end`: Genomic coordinates of exons.
#'   \item `width`: Length of each exon.
#'   \item `strand`: Strand direction (`+` or `-`).
#'   \item `gene_name`: Gene name.
#'   \item `transcript_id`: Transcript ID.
#'   \item `tx_len`: Cumulative transcript length up to this exon.
#'   \item `f_len`: Total length of the transcript.
#' }
#'
#' @importFrom dplyr filter group_by summarise arrange mutate left_join select
#'
#' @examples
#' \dontrun{
#' gtf_file <- "path/to/your.gtf"
#' transcripts <- prepareTransInfo_forGenome(gtf_file)
#' head(transcripts)
#' }
#'
#' @export
prepareTransInfo_forGenome <- function(gtf_file = NULL){
  # load gtf
  if(is.character(gtf_file)){
    if (requireNamespace("rtracklayer", quietly = TRUE)) {
      gtf <- rtracklayer::import.gff(gtf_file,format = "gtf") %>%
        data.frame()
    } else {
      warning("Package 'rtracklayer' is needed for this function to work.")
    }
  }else if("data.frame" %in% class(gtf_file)){
    gtf <- gtf_file
  }


  # total transcript lenth
  features_len <- gtf %>%
    dplyr::filter(type == "exon") %>%
    dplyr::group_by(gene_name,transcript_id) %>%
    dplyr::summarise(f_len = sum(width),.groups = "drop")

  # calculate transcript position block for positive strand
  exon_pos <- gtf %>%
    dplyr::filter(type == "exon" & strand == "+") %>%
    dplyr::select(seqnames,start,end,width,strand,gene_name,transcript_id) %>%
    dplyr::group_by(gene_name,transcript_id) %>%
    dplyr::arrange(start,end,.by_group = T) %>%
    dplyr::mutate(tx_len = cumsum(width),.after = "width") %>%
    dplyr::left_join(y = features_len,by = c("gene_name","transcript_id"))

  # calculate transcript position block for negtive strand
  exon_neg <- gtf %>%
    dplyr::filter(type == "exon" & strand == "-") %>%
    dplyr::select(seqnames,start,end,width,strand,gene_name,transcript_id) %>%
    dplyr::group_by(gene_name,transcript_id) %>%
    dplyr::arrange(dplyr::desc(start),dplyr::desc(end),.by_group = T) %>%
    dplyr::mutate(tx_len = cumsum(width),.after = "width") %>%
    dplyr::left_join(y = features_len,by = c("gene_name","transcript_id"))

  # combine + and - data
  dist_features <- rbind(exon_pos, exon_neg)

  return(dist_features)
}


# ==============================================================================
# function to calculate ribosome occupancy on genome
# ==============================================================================

#' Compute Read Occupancy on the Genome
#'
#' The `getOccupancyGenome()` function extracts read positions from a BAM file
#' for a given gene, converts them to transcriptome positions, and calculates
#' normalized read occupancy (RPM).
#'
#' @param bam_file Character. Path to the BAM file.
#' @param gene_name Character. The gene name for which read occupancy is extracted.
#' @param features A `data.frame` or `GRanges` object containing transcript annotation.
#' @param total_reads Numeric. The total number of mapped reads (for RPM normalization).
#' @param assignment_mode Character. Specifies the read assignment strategy: `"end5"` (5' end) or `"end3"` (3' end).
#' @param coordinate_to_trans Logical. Whether to convert genomic coordinates to transcript
#'        coordinates. Default is FALSE.
#'
#' @return A `data.frame` with read occupancy information:
#' \itemize{
#'   \item `rname`: Transcript and gene identifier (`"transcript_id|gene_name"`).
#'   \item `pos`: Transcriptomic read position.
#'   \item `count`: Number of reads at the position.
#'   \item `rpm`: Reads per million (RPM) normalized occupancy.
#' }
#'
#' @importFrom dplyr rename mutate group_by summarise filter select case_when
#' @importFrom fastplyr f_filter
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges findOverlaps
#' @importFrom Rsamtools scanBam ScanBamParam scanBamFlag
#' @importFrom parallel detectCores
#'
#' @examples
#' \dontrun{
#' bam_file <- "example.bam"
#' gene <- "TP53"
#' total_mapped_reads <- 1e7
#' features <- prepareTransInfo_forGenome(gtf_file = "example.gtf")
#'
#' occupancy_data <- getOccupancyGenome(
#'   bam_file = bam_file,
#'   gene_name = gene,
#'   features = features,
#'   total_reads = total_mapped_reads
#' )
#' head(occupancy_data)
#' }
#'
#' @export
getOccupancyGenome <- function(bam_file = NULL,
                               gene_name = NULL,
                               features = NULL,
                               total_reads = NULL,
                               assignment_mode = NULL,
                               coordinate_to_trans = FALSE){
  # filter genes
  query_region <- features %>%
    dplyr::rename(gene = gene_name) %>%
    fastplyr::f_filter(gene == gene_name) %>%
    GenomicRanges::GRanges()

  # get position
  bam_data <- Rsamtools::scanBam(file = bam_file,
                                 nThreads = parallel::detectCores(),
                                 param = Rsamtools::ScanBamParam(what = c("rname", "pos", "strand", "flag"),
                                                                 which = query_region,
                                                                 flag = Rsamtools::scanBamFlag(isUnmappedQuery = FALSE))
  )

  # list to data frame
  bminfo <- lapply(bam_data,FUN = data.frame) %>% do.call("rbind",.)

  if(assignment_mode == "end5"){
    # check strand to assign 5'end position
    bminfo <- bminfo %>%
      dplyr::mutate(pos = ifelse(strand == "-", pos + qwidth - 1, pos))
  }else{
    # check strand to assign 3'end position
    bminfo <- bminfo %>%
      dplyr::mutate(pos = ifelse(strand == "+", pos + qwidth - 1, pos))
  }

  bminfo <- bminfo %>%
    dplyr::group_by(rname,pos) %>%
    dplyr::summarise(count = dplyr::n(),.groups = "drop") %>%
    dplyr::filter(pos >= min(IRanges::start(query_region@ranges)) &
                    pos <= max(IRanges::end(query_region@ranges)))

  # ============================================================================
  # check coordinate transform
  if(coordinate_to_trans == TRUE){
    bminfo <- bminfo %>%
      dplyr::rename(seqnames = rname, start = pos) %>%
      dplyr::mutate(end = start,.after = "start")

    # to GRanges
    bminfo.rg <- GenomicRanges::GRanges(bminfo)

    # transcript annotation regions
    trans_rg <- GenomicRanges::GRanges(features)


    # genomic position to transcriptome position
    ov <- IRanges::findOverlaps(query = bminfo.rg,subject = trans_rg)

    # get overlap data
    if (requireNamespace("S4Vectors", quietly = TRUE)) {
      lo <- cbind(as.data.frame(bminfo.rg[S4Vectors::queryHits(ov)]),
                  as.data.frame(trans_rg[S4Vectors::subjectHits(ov)]))
    } else {
      warning("Package 'S4Vectors' is needed for this function to work.")
    }


    # make unique names
    names(lo) <- make.names(names(lo),unique = T)

    # calculate transcript position
    tgene <- gene_name
    lo <- lo %>%
      dplyr::mutate(pos = dplyr::case_when(strand.1 == "+" ~ tx_len - abs(end.1 - start),
                                           strand.1 == "-" ~ tx_len - abs(start - start.1)),
                    rname = paste(transcript_id,gene_name,sep = "|")) %>%
      dplyr::filter(gene_name == tgene)

    # output
    tdf <- lo %>% dplyr::select(rname,pos,count)
  }else{
    tdf <- bminfo
  }

  # rpm normalization
  tdf$rpm <- (tdf$count/total_reads)*10^6


  return(tdf)
}



# ==============================================================================
# function to calculate rna coverage on genome
# ==============================================================================

#' Compute RNA-seq Coverage on the Genome
#'
#' The `getCoverageGenome()` function extracts read coverage information from a BAM file
#' for a given gene, converts it to transcriptomic positions, and calculates normalized
#' coverage (RPM).
#'
#' @param bam_file Character. Path to the RNA-seq BAM file.
#' @param gene_name Character. The gene name for which transcript coverage is extracted.
#' @param features A `data.frame` or `GRanges` object containing transcript annotation.
#' @param total_reads Numeric. The total number of mapped reads (for RPM normalization).
#' @param coordinate_to_trans Logical. Whether to convert genomic coordinates to transcript
#'        coordinates. Default is FALSE.
#'
#' @return A `data.frame` with transcript coverage information:
#' \itemize{
#'   \item `rname`: Transcript and gene identifier (`"transcript_id|gene_name"`).
#'   \item `pos`: Transcriptomic position of the read coverage.
#'   \item `count`: Coverage (number of reads covering the position).
#'   \item `rpm`: Reads per million (RPM) normalized coverage.
#' }
#'
#' @importFrom dplyr rename mutate filter select case_when
#' @importFrom fastplyr f_filter
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges findOverlaps
#' @importFrom Rsamtools pileup BamFile PileupParam ScanBamParam
#'
#' @examples
#' \dontrun{
#' bam_file <- "example.bam"
#' gene <- "TP53"
#' total_mapped_reads <- 1e7
#' features <- prepareTransInfo_forGenome(gtf_file = "example.gtf")
#'
#' coverage_data <- getCoverageGenome(
#'   bam_file = bam_file,
#'   gene_name = gene,
#'   features = features,
#'   total_reads = total_mapped_reads
#' )
#' head(coverage_data)
#' }
#'
#' @export
getCoverageGenome <- function(bam_file = NULL,
                              gene_name = NULL,
                              features = NULL,
                              total_reads = NULL,
                              coordinate_to_trans = FALSE){
  # filter genes
  query_region <- features %>%
    dplyr::rename(gene = gene_name) %>%
    fastplyr::f_filter(gene == gene_name) %>%
    GenomicRanges::GRanges()

  # get coverage
  pileup_result <- Rsamtools::pileup(file = Rsamtools::BamFile(bam_file),
                                     pileupParam = Rsamtools::PileupParam(distinguish_nucleotides = FALSE,
                                                                          distinguish_strands = FALSE,
                                                                          max_depth = 10^6),
                                     scanBamParam = Rsamtools::ScanBamParam(which = query_region)) %>%
    dplyr::select(-which_label)

  colnames(pileup_result)[1] <- "rname"

  # ============================================================================
  # check coordinate transform
  if(coordinate_to_trans == TRUE){
    pileinfo <- pileup_result %>%
      dplyr::rename(seqnames = rname, start = pos) %>%
      dplyr::mutate(end = start,.after = "start")

    # to GRanges
    pileinfo.rg <- GenomicRanges::GRanges(pileinfo)

    # transcript annotation regions
    trans_rg <- GenomicRanges::GRanges(features)


    # genomic position to transcriptome position
    ov <- IRanges::findOverlaps(query = pileinfo.rg,subject = trans_rg)

    # get overlap data
    if (requireNamespace("S4Vectors", quietly = TRUE)) {
      lo <- cbind(as.data.frame(pileinfo.rg[S4Vectors::queryHits(ov)]),
                  as.data.frame(trans_rg[S4Vectors::subjectHits(ov)]))
    } else {
      warning("Package 'S4Vectors' is needed for this function to work.")
    }


    # make unique names
    names(lo) <- make.names(names(lo),unique = T)

    # calculate transcript position
    tgene <- gene_name
    lo <- lo %>%
      dplyr::mutate(pos = dplyr::case_when(strand.1 == "+" ~ tx_len - abs(end.1 - start),
                                           strand.1 == "-" ~ tx_len - abs(start - start.1)),
                    rname = paste(transcript_id,gene_name,sep = "|")) %>%
      dplyr::filter(gene_name == tgene)

    # output
    tdf <- lo %>%
      dplyr::select(rname,pos,count)
  }else{
    tdf <- pileup_result %>%
      dplyr::filter(pos >= min(IRanges::start(query_region@ranges)) &
                      pos <= max(IRanges::end(query_region@ranges)))
  }


  # rpm normalization
  tdf$rpm <- (tdf$count/total_reads)*10^6

  return(tdf)
}




# ==============================================================================
# function to visualize coverage and occupancy on genome
# ==============================================================================

#' @title Generic function for genome transcript plotting
#'
#' @description
#' This generic function serves as an interface for plotting genomic transcripts.
#' It allows specific methods to be applied depending on the class of the object provided.
#'
#' @param object An object of class \code{ribotrans}, or other classes with defined methods for genome transcript plotting.
#' @param ... Additional arguments passed to specific methods.
#'
#' @return A ggplot object or other relevant output depending on the method used.
#'
#' @export
#' @rdname genome_trans_plot
setGeneric("genome_trans_plot",function(object,...) standardGeneric("genome_trans_plot"))




#' @title Method for plotting genomic transcripts for ribotrans objects
#'
#' @description
#' This method generates genomic plots for ribosomal occupancy and RNA reads coverage
#' based on the specified parameters for objects of class \code{ribotrans}.
#'
#' @param object An object of class \code{ribotrans} containing the necessary data for plotting.
#' @param selected_id (optional) A specific transcript id for the gene of interest. Default is NULL.
#' @param type A character string indicating the type of data to plot. Choices are "ribo" for
#' ribosomal occupancy, "rna" for RNA coverage, or "ribo_rna" for both. Default is "ribo".
#' @param layer A character string indicating the plot layer type, either "col" for `geom_col` plots
#' or "line" for `geom_path` plots. Default is "col".
#' @param sample_order (optional) A character vector specifying the order in which samples should be
#' displayed in the plot. Default is NULL.
#' @param facet_layer A facet layer for ggplot, default is \code{ggplot2::facet_grid(sample~rname,switch = "y")}.
#' @param sep_mer_sample A boolean indicating whether to separate merged and individual samples. Default is FALSE.
#' @param new_signal_range A boolean indicating whether to show the new signal range on the plot. Default is FALSE.
#' @param collapse_structure A boolean indicating whether to collapse the structures of transcripts and exons into a single line in the plot. Default is FALSE.
#' @param range_x A numeric value specifying the position adjustment for the range label on the x-axis. Default is 0.9.
#' @param range_y A numeric value specifying the position adjustment for the range label on the y-axis. Default is 0.9.
#' @param range_size A numeric value specifying the font size for the range label. Default is 4.
#' @param range_digit An integer indicating the number of digits to round the range label. Default is 1.
#' @param scale_factor A scaling factor to adjust RNA coverage values for better visualization. Default is 1.
#' @param utr_width A numeric value to specify the width of the UTR segments in the plot. Default is 1.
#' @param cds_width A numeric value to specify the width of the CDS segments in the plot. Default is 3.
#' @param background_line A vector defining the color and linewidth for the background line structure. Default is \code{c("black", 0.25)}.
#' @param exon_line A vector defining the color and linewidth for the exon line structure. Default is \code{c("#003399", 1)}.
#' @param cds_line A vector defining the color and linewidth for the CDS line structure. Default is \code{c("#003399", 3)}.
#' @param structure_panel_height A numeric value specifying the height of the structure panel in the plot. Default is 0.2.
#'
#' @return A ggplot object displaying the genomic transcripts along with ribosome occupancy or RNA coverage data.
#'
#' @details The function uses ggplot2 for visualization and requires specific data to be pre-calculated (e.g., ribosome occupancy or RNA coverage). Users must ensure that the appropriate methods (such as \code{get_occupancy} or \code{get_coverage}) have been executed prior to plotting.
#'
#' @examples
#' \dontrun{
#' # Assuming 'ribotrans_obj' is a pre-defined ribotrans object
#' genome_trans_plot(
#'   object = ribotrans_obj,
#'   type = "ribo",
#'   layer = "col",
#'   sample_order = c("sample1", "sample2", "merged_sample"),
#'   collapse_structure = TRUE
#' )
#' }
#'
#' @import ggplot2
#' @import dplyr
#' @importFrom ggside scale_xsidey_continuous
#' @importFrom ggside ggside
#'
#' @export
#' @rdname genome_trans_plot
setMethod("genome_trans_plot",
          signature(object = "ribotrans"),
          function(object,
                   selected_id = NULL,
                   type = c("ribo","rna","ribo_rna"),
                   layer = c("col", "line"),
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
                player <- geom_col(aes(x = pos,y = smooth,fill = sample),width = 1)
              }

              col <- scale_fill_manual(values = rep("black",length(unique(pldf$sample))))
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
                player <- geom_col(aes(x = pos,y = smooth,fill = sample),width = 1)
              }

              col <- scale_fill_manual(values = rep("black",length(unique(pldf$sample))))
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
                col <- scale_color_manual(values = c("ribo" = "red", "rna" = "grey50"))
              }else{
                player <- geom_col(aes(x = pos,y = smooth,fill = exp),width = 1)
                col <- scale_fill_manual(values = c("ribo" = "red", "rna" = "grey50"))
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
              dplyr::filter(gene_name == gene & type == "transcript") %>%
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
              dplyr::filter(gene_name == gene & type == "exon") %>%
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
              dplyr::filter(gene_name == gene & type == "CDS") %>%
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
                tmp.df$sp <- factor(pldf$sp ,levels = c(sample_order,"merged sample"))
              }else{
                tmp.df$sample <- factor(pldf$sample ,levels = sample_order)
              }

            }

            # ==================================================================
            # plot
            # p <-
            ggplot(pldf) +
              player +
              range_label +
              # theme_bw() +
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
              col

          }
)

