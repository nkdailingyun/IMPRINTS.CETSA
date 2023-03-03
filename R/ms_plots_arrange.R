#' ms_plots_arrange
#'
#' Function to generate pdf files with multipanel compiled from plots objects
#'
#' @param plotlist a list of plots objects (now hardcode for the case of 2,3,4)
#' @param innerlayout a vector indicating the lower level panel layout for
#' each multi-panel plot, default value is NULL, would change according to
#' the number of supplied plots
#' @param mainlayout a vector indicating the main panel layout for multi-panel
#' plots per page, default value is c(4,1)
#' @param toplabel textual label at the top part of the page
#' @param leftlabel textual label at the left side of the page
#' @param bottomlabel textual label at the bottom part of the page
#' @param pdfheight a number indicate the height of pdf file, default value 12
#' @param pdfwidth a number indicate the width of pdf file, default value 12
#'
#'
#' @import RColorBrewer
#' @import ggplot2
#' @importFrom gridExtra grid.arrange arrangeGrob
#'
#' @export
#' @keywords internal
#'
#' @return a list of ggplot2 object
#' @examples \dontrun{
#'   ms_plots_arrange(MOLM)
#' }
#'
#'

ms_plots_arrange <- function(plotlist, innerlayout=NULL, mainlayout=c(4,1),
                             external=TRUE, toplabel="", leftlabel="", bottomlabel="",
                             pdfname="mixed_ggplotting.pdf",
                             pdfheight=12, pdfwidth=10, returnplots=FALSE) {

  if (external) { external_graphs(T) }
  outdir <- NULL
  for (i in seq_along(plotlist)) {
    assign(paste0("plot",i), plotlist[[i]])
  }
  plots <- list()
  names <- unique(names(plot1))
  if (length(plotlist)==2) {
    innerlayout = c(1,2)
    for (name in names) {
      plots[[name]] <- gridExtra::grid.arrange(plot1[[name]], plot2[[name]],
                                               ncol=innerlayout[2], nrow=innerlayout[1])
    }
  } else if (length(plotlist)==3) {
    innerlayout = c(1,3)
    for (name in names) {
      plots[[name]] <- gridExtra::grid.arrange(plot1[[name]], plot2[[name]], plot3[[name]],
                                               ncol=innerlayout[2], nrow=innerlayout[1])
    }
  } else if (length(plotlist)==4) {
    innerlayout = c(1,4)
    for (name in names) {
      plots[[name]] <- gridExtra::grid.arrange(plot1[[name]], plot2[[name]],
                                               plot3[[name]], plot4[[name]],
                                               ncol=innerlayout[2], nrow=innerlayout[1])
    }
  }
  # print(class(plots[[1]]))
  # print(class(plots))
  if (returnplots) { return(plots) }

  params <- list(nrow=mainlayout[1], ncol=mainlayout[2])
  n <- with(params, nrow*ncol)
  ## add one page if division is not complete
  pages <- length(plots) %/% n + as.logical(length(plots) %% n)
  groups <- split(seq_along(plots), gl(pages, n, length(plots)))

  pl <- lapply(names(groups), function(i){

    gridExtra::grid.arrange(
      do.call(gridExtra::arrangeGrob,
              c(plots[groups[[i]]], params, top=toplabel,
                left=leftlabel,bottom=bottomlabel)))
  })

  class(pl) <- c("arrangelist", "ggplot", class(pl))
  pdfname <- gsub("/", " ", pdfname)
  if (length(outdir)) {
    ggplot2::ggsave(file=paste0(outdir,"/",format(Sys.time(), "%y%m%d_%H%M_"), pdfname),
           pl, height=pdfheight, width=pdfwidth)
  } else {
    ggplot2::ggsave(file=paste0(format(Sys.time(), "%y%m%d_%H%M_"), pdfname), pl,
           height=pdfheight, width=pdfwidth)
  }

  if (external) { external_graphs(F) } # switch off the external graphs
  message("IMPRINTS-CETSA bar plot file generated successfully.")
}


