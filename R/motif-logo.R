# ==============================================================================
# logo plot
# ==============================================================================

#' Sequence Logo Plot for Biological Motifs
#'
#' This function generates a sequence logo visualization based on aligned biological
#' sequences using multiple methods, including Shannon information content, relative
#' frequency, EDLogo and log2 enrichment (foreground vs. background). It supports amino acid, DNA and RNA sequence types.
#'
#' @param foreground_seqs A character vector or Biostrings::XStringSet object of foreground sequences (e.g., motif-aligned peptides or DNA k-mers).
#' @param background_seqs A character vector or Biostrings::XStringSet of background sequences. Required when \code{method = "enrich"}.
#' @param rev_stack_order Logical. Reverse stacking order of letters. Default is \code{FALSE}.
#' @param type Character. Only applies to \code{method = "enrich"}. If \code{"merge"} (default),
#' enriched and depleted residues are shown in one logo with positive/negative y-axis. If \code{"sep"},
#' the two are split into two subplots (up/down).
#' @param col_scheme Character. Color scheme to be passed to \code{ggseqlogo}. Examples include \code{"chemistry"}, \code{"hydro"}, \code{"classic"}, etc.
#' @param seq_type Character. Type of sequence: \code{"aa"} (amino acids), \code{"dna"}, or \code{"rna"}.
#' @param method Character. One of \code{"bits"}, \code{"prob"}, \code{"EDLogo"}, or \code{"enrich"}:
#' \itemize{
#' \item \code{"bits"}: Shannon information content.
#' \item \code{"prob"}: Residue frequency logo.
#' \item \code{"EDLogo"}: Positional deviation using the Logolas package.
#' \item \code{"enrich"}: Log2 enrichment (foreground vs. background) at each position, similar to pLogo-style representation (not statistical).
#' }
#' @param return_data Logical. If \code{TRUE} and \code{method = "enrich"}, returns the numeric log2 enrichment matrix instead of a plot object.
#'
#' @return A ggplot object representing a sequence logo, or a numeric enrichment matrix if \code{return_data = TRUE} and \code{method = "enrich"}.
#'
#' @details
#' For \code{method = "enrich"}, the function computes:
#' \deqn{ \log_2 \left(\frac{P_\mathrm{fg}}{P_\mathrm{bg}} \right) }
#' where \eqn{P_\mathrm{fg}}, \eqn{P_\mathrm{bg}} are position-specific residue frequencies.
#' Residues with infinite log ratios are set to 0 to avoid plotting artifacts.
#'
#' Unlike pLogo, this method does not compute statistical significance (no binomial test or
#' Bonferroni adjustment), but provides a fast approximate visual signature.
#'
#' @examples
#' \dontrun{
#' library(Biostrings)
#' # Example: log2 enrichment between foreground and background
#' fg <- c("ARND", "ARNE", "ARNQ")
#' bg <- c("ARND", "GHIY", "PQRS", "ARND", "ARND")
#' logo_plot(foreground_seqs = fg, background_seqs = bg, method = "enrich", seq_type = "aa")
#'
#' # AA probabilistic logo
#' logo_plot(foreground_seqs = fg, method = "prob", seq_type = "aa")
#'
#' # DNA bits logo
#' dna_fg <- c("ATGC", "ATGT", "ATGA")
#' logo_plot(foreground_seqs = dna_fg, method = "bits", seq_type = "dna")
#' }
#'
#' @seealso \code{\link[ggseqlogo]{ggseqlogo}}, \code{\link[Logolas]{logomaker}}, \code{\link{to_position_matrix}}
#' @import ggplot2
#' @importFrom Biostrings AAStringSet DNAStringSet RNAStringSet
#' @export
logo_plot <- function(foreground_seqs = NULL,
                      background_seqs = NULL,
                      rev_stack_order = FALSE,
                      type = c("merge","sep"),
                      col_scheme = "chemistry2",
                      seq_type = c("aa", "dna", "rna"),
                      method = c("bits","prob","EDLogo","enrich"),
                      return_data = FALSE){
  method <- match.arg(method, choices = c("bits","prob","EDLogo","enrich"))
  seq_type <- match.arg(seq_type, choices = c("aa", "dna", "rna"))
  # ============================================================================
  # check method
  if(method %in% c("bits","prob")){
    if (requireNamespace("ggseqlogo", quietly = TRUE)) {
      logo <- ggseqlogo::ggseqlogo(foreground_seqs, method = method,
                                   rev_stack_order = rev_stack_order,
                                   col_scheme = col_scheme,
                                   seq_type = seq_type) +
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

      # check seq_type
      if(seq_type == "aa"){
        fg_seq <- Biostrings::AAStringSet(foreground_seqs)
        bg_seq <- Biostrings::AAStringSet(background_seqs)
      }else if(seq_type == "dna"){
        fg_seq <- Biostrings::DNAStringSet(foreground_seqs)
        bg_seq <- Biostrings::DNAStringSet(background_seqs)
      }else if(seq_type == "rna"){
        fg_seq <- Biostrings::RNAStringSet(foreground_seqs)
        bg_seq <- Biostrings::RNAStringSet(background_seqs)
      }

      fg_matrix <- to_position_matrix(fg_seq, seq_type = seq_type)
      bg_matrix <- to_position_matrix(bg_seq, seq_type = seq_type)

      fg_matrix_freq <- sweep(fg_matrix, 2, colSums(fg_matrix), "/")
      bg_matrix_freq <- sweep(bg_matrix, 2, colSums(bg_matrix), "/")

      log2ratio <- log2(fg_matrix_freq/bg_matrix_freq)

      log2ratio[is.infinite(log2ratio)] <- 0

      # ========================================================================
      # plot
      if(type == "merge"){
        if (requireNamespace("ggseqlogo", quietly = TRUE)) {
          logo <- ggseqlogo::ggseqlogo(log2ratio,
                                       seq_type = seq_type,
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
                                 seq_type = seq_type,
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
                                 seq_type = seq_type,
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

#' Statistical Sequence Logo Plot (pLogo-style)
#'
#' This function generates a statistical sequence logo (pLogo-style) that visualizes overrepresented and
#' underrepresented residues based on a binomial test. It compares position-specific counts
#' in a foreground and background dataset using a log-odds scoring approach. A Bonferroni-corrected
#' significance threshold (P < 0.05) is used to display a reference line.
#'
#' @param foreground_seqs A character vector or Biostrings::XStringSet of aligned foreground sequences
#' (e.g. amino acid motif windows).
#' @param background_seqs A character vector or Biostrings::XStringSet of aligned background sequences. Required.
#' @param rev_stack_order Logical. Reverse the stacking order of letters (default: TRUE).
#' @param seq_type Character. Sequence type, one of \code{"aa"} (amino acid), \code{"dna"}, or \code{"rna"}.
#' @param type Character. Display style. If \code{"merge"}, plot shows both enriched/depleted residues
#' in one panel. If \code{"sep"}, returns a top and bottom panel separately (like canonical pLogo).
#' @param col_scheme Character. Color scheme to use, passed to \code{ggseqlogo}. Default: \code{"chemistry2"}.
#' @param npcx Numeric `[0–1]`. X position (NPC units) for annotation text.
#' @param npcy Numeric `[0–1]`. Y position (NPC units) for annotation text.
#' @param size Numeric. Text size of the annotation.
#' @param return_data Logical. If \code{TRUE}, returns the internally computed log-odds matrix instead of a plot.
#'
#' @return A ggplot object showing a sequence logo with log-probability-scaled residue heights, or a numeric
#' matrix of log-odds ratios if \code{return_data=TRUE}.
#'
#' @details
#' This implementation follows the principles of the pLogo algorithm (O’Shea et al, Nature Methods, 2013):
#' residues are scaled according to the log odds of their binomial probability of enrichment or depletion
#' given the background. Significant thresholds are marked with dashed lines, using Bonferroni-corrected tests.
#'
#' A red dashed horizontal line is shown at:
#' \deqn{\pm \log_{10}\left(\frac{\alpha'}{1-\alpha'}\right)}
#' where \eqn{\alpha'} is the corrected significance level, accounting for multiple testing:
#' \deqn{\alpha' = \frac{0.05}{N_{positions} \times N_{residues} - N_{zero\ fg}}}
#'
#' @references
#' O’Shea JP, Chou MF, Quader SA, Ryan JK, Church GM, Schwartz D. pLogo: a probabilistic
#' approach to visualizing sequence motifs. Nat Methods. 2013 Dec;10(12):1211–1212. \doi{10.1038/nmeth.2646}
#'
#'
#' @examples
#' \dontrun{
#' library(Biostrings)
#' fg <- c("ARND", "ARNE", "ARNQ", "ARNA")
#' bg <- rep(c("ARND", "GHIL", "PQRS"), each = 5)
#' plogo_plot(fg, bg, seq_type = "aa", type = "merge")
#' }
#'
#' @import ggplot2
#' @importFrom Biostrings AAStringSet DNAStringSet RNAStringSet
#' @export
plogo_plot <- function(foreground_seqs = NULL,
                       background_seqs = NULL,
                       rev_stack_order = TRUE,
                       seq_type = c("aa", "dna", "rna"),
                       type = c("merge","sep"),
                       col_scheme = "chemistry2",
                       npcx = 0.98, npcy = 0.98, size = 3.5,
                       return_data = FALSE){
  seq_type <- match.arg(seq_type, choices = c("aa", "dna", "rna"))
  type <- match.arg(type,choices = c("merge","sep"))
  # ============================================================================
  # enrichment logo
  # ==========================================================================

  # check seq_type
  if(seq_type == "aa"){
    fg_seq <- Biostrings::AAStringSet(foreground_seqs)
    bg_seq <- Biostrings::AAStringSet(background_seqs)
  }else if(seq_type == "dna"){
    fg_seq <- Biostrings::DNAStringSet(foreground_seqs)
    bg_seq <- Biostrings::DNAStringSet(background_seqs)
  }else if(seq_type == "rna"){
    fg_seq <- Biostrings::RNAStringSet(foreground_seqs)
    bg_seq <- Biostrings::RNAStringSet(background_seqs)
  }

  fg_matrix <- to_position_matrix(fg_seq, seq_type = seq_type)
  bg_matrix <- to_position_matrix(bg_seq, seq_type = seq_type)

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
                           seq_type = seq_type,
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
                             seq_type = seq_type,
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
                             seq_type = seq_type,
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
