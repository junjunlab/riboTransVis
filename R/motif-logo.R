# ==============================================================================
# logo plot
# ==============================================================================

#' Plot Sequence Logos with Multiple Visualization Methods
#'
#' This function creates sequence logo plots using various visualization approaches,
#' including Shannon bits, probabilities, EDLogo, and enrichment-based log2-ratio logos.
#'
#' @param foreground_seqs A character vector of foreground sequences (aligned peptide motifs) or a matrix. This input is used across all plotting modes.
#' @param background_seqs Optional character vector of background sequences. Only used when \code{method = "enrich"} to compute enrichment over a background model.
#' @param rev_stack_order Logical. Whether to reverse the stacking order of residues within each logo column. Default is \code{FALSE}.
#' @param type Character. Either \code{"merge"} (default) to plot a unified enrichment logo or \code{"sep"} to plot enriched and depleted residues in separate panels. Only applies when \code{method = "enrich"}.
#' @param col_scheme Character. Amino acid color scheme used for rendering logos. Passed to \code{ggseqlogo()}.
#' @param method Character. Visualization method. Can be one of:
#' \itemize{
#'   \item \code{"bits"}: Shannon bit-based sequence logo (information content, uses \pkg{ggseqlogo}).
#'   \item \code{"prob"}: Residue probability-based logo (uses \pkg{ggseqlogo}).
#'   \item \code{"EDLogo"}: Uses EDLogo rendering from package \pkg{Logolas}.
#'   \item \code{"enrich"}: Enrichment logo showing log2 ratio of foreground vs background residue frequencies.
#' }
#' @param return_data Logical. If \code{TRUE}, returns the log2 enrichment matrix instead of a plot (only valid when \code{method = "enrich"}). Defaults to \code{FALSE}.
#'
#' @return A \code{ggplot} object (or a \code{cowplot} object when \code{type = "sep"}) representing the motif logo visualization.
#'
#' @details
#' \code{logo_plot()} allows for visualizing motifs using traditional or comparative methods:
#' \itemize{
#'   \item For \code{"bits"} and \code{"prob"}, only foreground data is used.
#'   \item For \code{"EDLogo"}, a matrix or data frame is expected and will be passed to \code{Logolas::logomaker()}.
#'   \item For \code{"enrich"}, foreground and background sequences are processed into position-specific frequency matrices, followed by computation of a log2-ratio enrichment matrix. This matrix is used to construct the pLogo-style enriched/depleted visualization.
#'   \item Negative log2 values indicate depletion; positive values indicate enrichment.
#' }
#' In \code{"enrich"} mode with \code{type = "sep"}, the plot is split into "enriched" (positive) and "depleted" (negative) panels using \pkg{cowplot}.
#'
#' @note The \code{"enrich"} method requires background sequences and assumes sequences are aligned and of the same length. \code{Biostrings} is used internally for consensus matrix generation.
#'
#' @importFrom Biostrings AAStringSet consensusMatrix
#' @importFrom ggplot2 theme_bw theme element_blank element_text element_line ylab geom_hline
#'
#' @examples
#' \dontrun{
#' ## Bits-based logo
#' fg <- c("AKTGRRKS", "AKSGRRKS", "AKTRRRRS")
#' logo_plot(foreground_seqs = fg, method = "bits")
#'
#' ## Enrichment logo (merge)
#' bg <- c("AKAAAKTS", "AKVVVKTS", "AKGGGKTS", "AKSSSKTS")
#' logo_plot(foreground_seqs = fg, background_seqs = bg, method = "enrich", type = "merge")
#'
#' ## Enrichment logo (split panels)
#' logo_plot(foreground_seqs = fg, background_seqs = bg, method = "enrich", type = "sep")
#' }
#'
#' @export
logo_plot <- function(foreground_seqs = NULL,
                      background_seqs = NULL,
                      rev_stack_order = FALSE,
                      type = c("merge","sep"),
                      col_scheme = "chemistry2",
                      method = c("bits","prob","EDLogo","enrich"),
                      return_data = FALSE){
  method <- match.arg(method, c("bits","prob","EDLogo","enrich"))
  # ============================================================================
  # check method
  if(method %in% c("bits","prob")){
    if (requireNamespace("ggseqlogo", quietly = TRUE)) {
      logo <- ggseqlogo::ggseqlogo(foreground_seqs, method = method,
                                   rev_stack_order = rev_stack_order,
                                   col_scheme = col_scheme) +
        theme_bw() +
        theme(panel.grid = element_blank(),
              axis.text = element_text(color = "black"))
    } else {
      warning("Package 'ggseqlogo' is needed for this function to work.")
    }

  }else if(method == "EDLogo"){
    if (requireNamespace("Logolas", quietly = TRUE)) {
      logo <- Logolas::logomaker(foreground_seqs, type = "EDLogo")
    } else {
      warning("Package 'Logolas' is needed for this function to work.")
    }
  }else{
    # ==========================================================================
    # enrichment logo
    # ==========================================================================
    if(!is.null(background_seqs)){
      to_position_matrix <- function(seqs){
        AA_STANDARD <- c("A","R","N","D","C","E","Q","G","H","I","L",
                         "K","M","F","P","S","T","W","Y","V")

        mat <- Biostrings::consensusMatrix(AAStringSet(seqs), as.prob = FALSE)
        mat <- mat[intersect(AA_STANDARD, rownames(mat)), ]
        return(mat)
      }

      fg_matrix <- to_position_matrix(Biostrings::AAStringSet(foreground_seqs))
      bg_matrix <- to_position_matrix(Biostrings::AAStringSet(background_seqs))

      fg_matrix_freq <- sweep(fg_matrix, 2, colSums(fg_matrix), "/")
      bg_matrix_freq <- sweep(bg_matrix, 2, colSums(bg_matrix), "/")

      log2ratio <- log2(fg_matrix_freq/bg_matrix_freq)

      log2ratio[is.infinite(log2ratio)] <- 0

      # ========================================================================
      # plot
      if(type == "merge"){
        if (requireNamespace("ggseqlogo", quietly = TRUE)) {
          logo <- ggseqlogo::ggseqlogo(log2ratio,
                                       seq_type == "aa",
                                       method = "custom",
                                       rev_stack_order = rev_stack_order,
                                       col_scheme = col_scheme) +
            geom_hline(yintercept = 0,lty = "solid", color = "red", linewidth = 1.25) +
            theme_bw() +
            theme(panel.grid = element_blank(),
                  axis.text = element_text(color = "black")) +
            ylab("Log2 enrichment ratio of probability")
        } else {
          warning("Package 'ggseqlogo' is needed for this function to work.")
        }
      }else{
        if (requireNamespace(c("ggseqlogo","cowplot","ggpp"), quietly = TRUE)){
          # filter positive and negtive residues
          log_score_mat_pos <- log2ratio
          log_score_mat_pos[log_score_mat_pos < 0] <- 0

          log_score_mat_neg <- log2ratio
          log_score_mat_neg[log_score_mat_neg > 0] <- 0

          # plot

          up <-
            ggseqlogo::ggseqlogo(log_score_mat_pos, method= "custom",
                                 rev_stack_order = rev_stack_order, col_scheme = col_scheme) +
            ggtitle("Log2 ratio of the frequency") +
            theme(legend.position = "none",
                  plot.title = element_text(hjust = 0.5),
                  axis.line.y = element_line(),
                  axis.ticks.y = element_line(),
                  axis.title.x = element_blank()) +
            ylab("Enriched")

          down <-
            ggseqlogo::ggseqlogo(log_score_mat_neg, method= "custom",
                                 rev_stack_order = rev_stack_order, col_scheme = col_scheme) +
            theme(axis.text.x = element_blank(),
                  axis.line.y = element_line(),
                  axis.ticks.y = element_line(),
                  plot.margin = margin(t = 0,unit = "cm")) +
            ylab("Depletted")

          # combine
          # logo <- cowplot::plot_grid(plotlist = list(up,down),align = "hv",ncol = 1)
          logo <- up/down
        }else{
          warning("Package 'ggseqlogo', 'cowplot' and 'ggpp' is needed for this function to work.")
        }


      }

    }
  }

  # return
  if(return_data == FALSE){
    return(logo)
  }else{
    return(log2ratio)
  }

}





# ==============================================================================
# define plogo plot
# ==============================================================================

#' Plot a pLogo-style Logo Based on Foreground and Background Sequences
#'
#' This function visualizes statistically enriched or depleted amino acid positions
#' from aligned sequences. It uses binomial probabilities to calculate log-odds scores
#' and generates a pLogo-like sequence logo using \code{ggseqlogo}. Both enriched and
#' depleted residues are shown, and a Bonferroni-corrected significance bar is added.
#'
#' @param foreground_seqs A character vector of amino acid sequences representing the foreground dataset.
#' @param background_seqs A character vector of amino acid sequences representing the background dataset.
#' @param rev_stack_order Logical. Should the residue stacks be ordered with the most significant on top? Default is TRUE.
#' @param type Character. Choose between \code{"merge"} (combined enriched and depleted in one plot) or \code{"sep"} (two-panel for enriched/depleted). Default is "merge".
#' @param col_scheme Residue color scheme passed to \code{ggseqlogo}; default is \code{"chemistry2"}.
#' @param npcx Numeric. Relative x-position for annotation text via \code{geom_text_npc()}; default is 0.98.
#' @param npcy Numeric. Relative y-position for annotation text; default is 0.98.
#' @param size Numeric. Font size for annotation text. Default is 3.5.
#' @param return_data Logical; if \code{TRUE}, returns the log-odds matrix instead of the plot. Default: \code{FALSE}.
#'
#' @return A ggplot object displaying the pLogo-style sequence logo, or a numeric matrix of log-odds enrichment scores if \code{return_data = TRUE}.
#'
#'
#' @details
#' Internally, this function performs several steps:
#' \itemize{
#'   \item Parses and counts each amino acid by position for both foreground and background sequences.
#'   \item Builds position matrices and computes frequency matrix for the background.
#'   \item Calculates residue log-odds using a binomial over/under probability ratio, via \code{get_log_ratio_mat()}.
#'   \item Adds red dashed lines to show Bonferroni-corrected p-value threshold (typically P < 0.05).
#'   \item Optionally splits plots into enriched and depleted panels (\code{type = "sep"}).
#' }
#'
#' Visualization is constructed using \code{ggseqlogo} and optionally arranged using the \code{patchwork} package.
#'
#' @importFrom Biostrings AAStringSet consensusMatrix
#' @importFrom ggplot2 geom_hline aes ylab theme element_text element_line
#'
#' @examples
#' \dontrun{
#' # Example: Generate logo from short protein sequences
#' fg <- c("ARGGGKTSV", "ARRRGKTSV", "ARSGGKTSV")
#' bg <- c("AKGGGKTSV", "ALGGGKTSV", "AKGGGKTSV", "ARAGGKTSV")
#'
#' plogo_plot(
#'   foreground_seqs = fg,
#'   background_seqs = bg,
#'   type = "merge",
#'   col_scheme = "chemistry2"
#' )
#' }
#'
#' @export
plogo_plot <- function(foreground_seqs = NULL,
                       background_seqs = NULL,
                       rev_stack_order = TRUE,
                       type = c("merge","sep"),
                       col_scheme = "chemistry2",
                       npcx = 0.98, npcy = 0.98, size = 3.5,
                       return_data = FALSE){
  type <- match.arg(type,choices = c("merge","sep"))
  # ============================================================================
  # enrichment logo
  # ==========================================================================

  to_position_matrix <- function(seqs){
    AA_STANDARD <- c("A","R","N","D","C","E","Q","G","H","I","L",
                     "K","M","F","P","S","T","W","Y","V")

    mat <- Biostrings::consensusMatrix(AAStringSet(seqs), as.prob = FALSE)
    mat <- mat[intersect(AA_STANDARD, rownames(mat)), ]
    return(mat)
  }

  fg_matrix <- to_position_matrix(Biostrings::AAStringSet(foreground_seqs))
  bg_matrix <- to_position_matrix(Biostrings::AAStringSet(background_seqs))

  # fg_matrix_freq <- sweep(fg_matrix, 2, colSums(fg_matrix), "/")
  bg_matrix_freq <- sweep(bg_matrix, 2, colSums(bg_matrix), "/")

  # get logp matrix
  log_score_mat <- get_log_ratio_mat(fg = fg_matrix,
                                     bg = bg_matrix,
                                     bg_freq = bg_matrix_freq)

  # get max height
  lapply(1:ncol(log_score_mat),function(x){
    tmp <- log_score_mat[,x]
    pos <- sum(tmp[tmp > 0])
    neg <- abs(sum(tmp[tmp < 0]))
    return(c(pos, neg))
  }) %>% unlist() %>% max() -> max.pos

  # assign max.pos
  log_score_mat[which(log_score_mat == 10)] <- max.pos

  # ============================================================================
  # check NULL amino acids in matrix
  null.n <- as.numeric(table(fg_matrix == 0)["TRUE"])

  # statistical-significance bar
  alph <- 0.05 / (ncol(fg_matrix)*nrow(fg_matrix) - null.n)
  ssbp <- log10(alph / (1 - alph))

  anno1 <- paste("(+|-)",round(abs(ssbp),digits = 2)," P < 0.05 n(fg) =",length(foreground_seqs))
  anno2 <- paste("fixed position n(bg) =",length(background_seqs))

  anno <- paste(anno1,anno2,sep = "\n")
  # ============================================================================
  # plot
  if (requireNamespace(c("ggseqlogo","cowplot","ggpp"), quietly = TRUE)) {
    if(type == "merge"){
      logo <-
      ggseqlogo::ggseqlogo(log_score_mat,method= "custom",
                rev_stack_order = rev_stack_order, col_scheme = col_scheme) +
        geom_hline(yintercept = 0,lty = "solid", color = "black", linewidth = 1.25) +
        geom_hline(yintercept = c(-ssbp,ssbp),lty = "dashed", color = "red", linewidth = 1) +
        ggpp::geom_text_npc(aes(label = anno,npcx = npcx,npcy = npcy),size = size) +
        theme(legend.position = "none",
              plot.title = element_text(hjust = 0.5),
              axis.line.y = element_line(),
              axis.ticks.y = element_line(),
              axis.title.x = element_blank()) +
        ylab("Log odds of the binomial probability")
    }else{
      # filter positive and negtive residues
      log_score_mat_pos <- log_score_mat
      log_score_mat_pos[log_score_mat_pos < 0] <- 0

      log_score_mat_neg <- log_score_mat
      log_score_mat_neg[log_score_mat_neg > 0] <- 0

      # plot

      up <-
        ggseqlogo::ggseqlogo(log_score_mat_pos, method= "custom",
                  rev_stack_order = rev_stack_order, col_scheme = col_scheme) +
        geom_hline(yintercept = -ssbp,lty = "dashed", color = "red", linewidth = 1) +
        ggpp::geom_text_npc(aes(label = anno,npcx = npcx,npcy = npcy),size = size) +
        ggtitle("Log odds of the binomial probability") +
        theme(legend.position = "none",
              plot.title = element_text(hjust = 0.5),
              axis.line.y = element_line(),
              axis.ticks.y = element_line(),
              axis.title.x = element_blank()) +
        ylab("Enriched")

      down <-
        ggseqlogo::ggseqlogo(log_score_mat_neg, method= "custom",
                  rev_stack_order = rev_stack_order, col_scheme = col_scheme) +
        geom_hline(yintercept = ssbp,lty = "dashed", color = "red", linewidth = 1) +
        theme(axis.text.x = element_blank(),
              axis.line.y = element_line(),
              axis.ticks.y = element_line(),
              plot.margin = margin(t = 0,unit = "cm")) +
        ylab("Depletted")

      # combine
      # logo <- cowplot::plot_grid(plotlist = list(up,down),align = "hv",ncol = 1)
      logo <- up/down

    }
  } else {
    warning("Package 'ggseqlogo', 'cowplot' and 'ggpp' are needed for this function to work.")
  }

  # return
  if(return_data == FALSE){
    return(logo)
  }else{
    return(log_score_mat)
  }
}
