#' imprints_barplotting
#'
#' Function to generate pdf files with multipanel bar plots for IMPRINTS-CETSA data
#'
#' @param data dataset after imprints_caldiff() to plot
#' @param treatmentlevel a vector of treatment labels, such as c("DMSO","TNFa","AT26533")
#' the order determines the arrangement, so in this case DMSO group would be the first group
#' @param setlevel a vector of set information if any, such as c("M13","M16")
#' @param corrtable a correlation table, could be useful for customized ranking
#' @param plotseq a vector of plots arrangement sequence (in composite ID)
#' @param printBothName a logical to tell if you want to print the both protein names on the plot
#' @param printGeneName a logical to tell if you want to print the gene names on the plot
#' @param pfdatabase a logical for using P.falciparum database or not
#' @param witherrorbar a logical to print or not the error bar on the plot
#' @param colorpanel a vector of customizable color scheme provided by the user, default set
#' c("gray","blue","orange")
#' @param usegradient whether the barplot should be draw in color gradient format
#' @param colorgradient the color scheme of gradient applied, default value c("#4575B4","ivory", "#D73027")
#' @param linegraph whether to plot the graph in a line graph format, default set to FALSE
#' @param log2scale whether the yscales should be in log2 scale, default set to TRUE
#' @param ratio aspect ratio of the plot, default set to 0.6
#' @param layout a vector indicating the panel layout for multi-panel plots
#' per page, default value is c(2,3) for set containing data, otherwise c(4,3)
#' @param ret_plot a logical to tell if you want to return the plots object
#' @param save_pdf a logical to tell if you want to save plots in a pdf file
#' @param toplabel textual label at the top part of the page
#' @param leftlabel textual label at the left side of the page
#' @param bottomlabel textual label at the bottom part of the page
#' @param pdfheight a number indicate the height of pdf file, default value 12
#' @param pdfwidth a number indicate the width of pdf file, default value 12
#'
#'
#' @importFrom tidyr expand gather separate unite
#' @importFrom dplyr group_by mutate rowwise summarise ungroup
#' @import ggplot2
#'
#' @export
#'
#' @return a list of ggplot2 object
#' @examples \dontrun{
#'   imprints_barplotting(MOLM, treatmentlevel=c("DMSO","TNFa","AT26533"), setlevel=c("M13","M16"))
#' }
#'
#'

imprints_barplotting <- function(data, treatmentlevel=NULL, setlevel=NULL, corrtable=NULL, plotseq=NULL,
                              printBothName=TRUE, printGeneName=FALSE, pfdatabase=FALSE,
                              witherrorbar=TRUE, colorpanel=c("gray","blue","orange"),
                              usegradient=FALSE, colorgradient=c("#4575B4","ivory", "#D73027"),
                              linegraph=FALSE, log2scale=TRUE, ratio=0.6, layout=NULL,
                              ret_plot=FALSE, save_pdf = TRUE, external=TRUE,
                              toplabel="IMPRINTS-CETSA bar plotting", leftlabel="", bottomlabel="",
                              pdfname="bar_ggplotting.pdf", pdfheight=12, pdfwidth=12) {

  # legenddata is any dataset containing the full levels of conditions, same as data
  dataname <- deparse(substitute(data))
  outdir <- ms_directory(data, dataname)$outdir
  data <- ms_directory(data, dataname)$data

  nrowdata <- nrow(data)
  if ( nrowdata==0 ) {
    message("Make sure there are more than one experimental condition in dataset.")
    stop("Otherwise specify remsinglecondprot==FALSE !")
  }

  # to concatenate id and description
  data <- combineProteinGeneName(data, printBothName, printGeneName, pfdatabase)

  data1 <- tidyr::gather(data[ ,!(names(data) %in% c("sumUniPeps","sumPSMs","countNum"))], condition, reading, -id)
  if (!log2scale) {
    data1 <- dplyr::mutate(data1, reading=2^reading)
  }
  a <- data1$condition[1]
  if (length(unlist(strsplit(a, "_")))==4) {
    withset <- TRUE
    data1 <- tidyr::separate(data1, condition, into=c("set","temperature","replicate","treatment"), sep="_")
    temperature <- sort(unique(data1$temperature))
    cdata <- plyr::ddply(data1, c("id", "set", "temperature", "treatment"), summarise,
                         N    = length(na.omit(reading)),
                         mean = mean(reading, na.rm=T),
                         sd   = sd(reading, na.rm=T),
                         se   = sd / sqrt(N)
    )
    if (length(layout)==0) {layout <- c(2,3)}
  } else if (length(unlist(strsplit(a, "_")))==3) {
    withset <- FALSE
    data1 <- tidyr::separate(data1, condition, into=c("temperature","replicate","treatment"), sep="_")
    temperature <- sort(unique(data1$temperature))
    cdata <- plyr::ddply(data1, c("id", "temperature", "treatment"), summarise,
                         N    = length(na.omit(reading)),
                         mean = mean(reading, na.rm=T),
                         sd   = sd(reading, na.rm=T),
                         se   = sd / sqrt(N)
    )
    if (length(layout)==0) {layout <- c(4,3)}
  } else {
    stop("make sure the namings of the columns of the dasaset are correct.")
  }

  cdata <- cdata %>% dplyr::rowwise() %>% dplyr::mutate(condition=paste(temperature, treatment, sep="_"))
  if (withset) { cdata$set <- factor(as.character(cdata$set), levels=setlevel) }
  if (class(corrtable)!="NULL") {
    corrtable <- corrtable[order(corrtable$correlation,decreasing=T), ]
    if (printBothName & !pfdatabase) {
      corrtable <- ms_composite_ID_Gene_Protein(corrtable,pfdatabase)
    } else if (printGeneName & !pfdatabase) {
      corrtable <- ms_composite_ID_Gene(corrtable,pfdatabase)
    } else {
      corrtable <- ms_composite_ID_Protein(corrtable,pfdatabase)
    }
    cdata$id <- factor(cdata$id, levels=corrtable$id)
  }
  if (length(plotseq)) {
    cdata$id <- factor(cdata$id, levels=plotseq)
  } else {
    cdata$id <- factor(cdata$id)
  }
  cdata$treatment <- factor(as.character(cdata$treatment), levels=treatmentlevel)
  cdata$condition <- factor(as.character(cdata$condition),
                            levels=apply(expand.grid(temperature,treatmentlevel), 1, paste, collapse="_"))

  message("Generating fitted plot file, pls wait.")
  # return(cdata)
  # print(head(cdata))

  if (external) { external_graphs(T) }

  plots <- plyr::dlply(cdata, plyr::.(id), .fun=ms_innerbarplotting, withset=withset)
  if (ret_plot) { return(plots) }
  params <- list(nrow=layout[1], ncol=layout[2])
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
  # print(class(pl)) #list
  class(pl) <- c("arrangelist", "ggplot", class(pl))
  pdfname <- gsub("/", " ", pdfname)
  if (length(outdir)) {
    ggplot2::ggsave(file=paste0(outdir,"/",format(Sys.time(), "%y%m%d_%H%M_"), dataname, "_", pdfname),
           pl, height=pdfheight, width=pdfwidth)
  } else {
    ggplot2::ggsave(file=paste0(format(Sys.time(), "%y%m%d_%H%M_"), dataname, "_", pdfname), pl,
           height=pdfheight, width=pdfwidth)
  }

  if (external) { external_graphs(F) } # switch off the external graphs
  message("IMPRINTS-CETSA bar plot file generated successfully.")
}
