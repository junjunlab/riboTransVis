#' Pipe operator
#'
#' See \code{magrittr::\link[magrittr:pipe]{\%>\%}} for details.
#'
#' @name %>%
#' @rdname pipe
#' @keywords internal
#' @export
#' @importFrom magrittr %>%
#' @usage lhs \%>\% rhs
#' @param lhs A value or the magrittr placeholder.
#' @param rhs A function call using the magrittr semantics.
#' @return The result of calling `rhs(lhs)`.
NULL




compact <- getFromNamespace("compact", "ggpp")
new_data_frame <- getFromNamespace("new_data_frame", "ggpp")
compact <- getFromNamespace("as_tibble", "tibble")

# globalVariables(c('geom', 'x', 'y', 'xmin', 'xmax', 'ymin', 'ymax', 'xend', 'yend', 'npcx',
#                   'npcy', 'label', '...', 'na.rm'))

annotate <-
  function (geom, x = NULL, y = NULL, xmin = NULL, xmax = NULL,
            ymin = NULL, ymax = NULL, xend = NULL, yend = NULL,
            npcx = NULL, npcy = NULL, label = NULL, ...,
            na.rm = FALSE)
  {
    if (inherits(label, what = c("data.frame", "gg", "grob"))) {
      label <- list(label)
    }
    position <- compact(list(x = x,
                             xmin = xmin,
                             xmax = xmax,
                             xend = xend,
                             y = y,
                             ymin = ymin,
                             ymax = ymax,
                             yend = yend,
                             npcx = npcx,
                             npcy = npcy,
                             label = label))
    aesthetics <- c(position, list(...))
    lengths <- vapply(aesthetics, length, integer(1))
    n <- unique(lengths)
    if (length(n) > 1L) {
      n <- setdiff(n, 1L)
    }
    if (length(n) > 1L) {
      bad <- lengths != 1L
      details <- paste(names(aesthetics)[bad], " (", lengths[bad],
                       ")", sep = "", collapse = ", ")
      rlang::abort(glue::glue("Unequal parameter lengths: {details}"))
    }
    data <- new_data_frame(position, n = n)
    ggplot2::layer(geom = geom, params = list(na.rm = na.rm, ...),
                   stat = StatIdentity,
                   position = PositionIdentity, data = data,
                   mapping = aes_all(names(data)),
                   inherit.aes = FALSE, show.legend = FALSE)
  }
