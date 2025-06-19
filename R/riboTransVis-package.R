#' @importFrom utils packageDescription
.onAttach <- function(libname, pkgname) {
  pkgVersion <- packageDescription(pkgname, fields = "Version")

  packageStartupMessage("Welcome to use riboTransVis for Ribo-seq analysis.")
  packageStartupMessage(paste("The version of riboTransVis:",
                              pkgVersion,
                              "\nAny advice or suggestions please contact with me: 3219030654@stu.cpu.edu.cn.",
                              sep = " "))
}
