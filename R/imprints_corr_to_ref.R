#' imprints_corr_in_complex
#'
#' Function to calculate the profile correlation between the subunit proteins of
#' the protein complexes found in IMPRINTS dataset
#'
#' @param data dataset after imprints_complex_mapping()
#' @param nread number of reading channels or sample treatments, default value 6
#' @param goodcorrcutoff the threshold for a good correlation, default value is 0.5
#' @param vgoodcorrcutoff the threshold for a very good correlation, default value is 0.9
#' @param removeredundancy scrutinize the complex input, to remove the redundancy
#' according to the measured subunits, default set to FALSE
#' @param similaritythreshold the threshold for similarity, numeric number between 0 to 1,
#' default use 1.0 to only remove the complete overlap ones
#' @param complexID a vector of complexID to retrieve from the data, if specified
#' @param doplotting whether to generate the correlation plot, default set to TRUE
#' @param layout a vector indicating the panel layout for multi-panel plots per page
#'
#'
#' @importFrom plyr . ddply
#' @importFrom ggplot2 ggsave
#' @importFrom GGally ggcorr
#' @export
#' @return a dataframe
#' @examples \dontrun{
#'   MOLM <- imprints_corr_in_complex(MOLM, nread=6)
#' }
#'
#'

imprints_corr_in_complex <- function(data=NULL, nread=6, goodcorrcutoff=0.5, vgoodcorrcutoff=0.9,
                                  removeredundancy=FALSE, similaritythreshold=1.0, complexID=NULL,
                                  doplotting=TRUE, layout=c(3,2), external=TRUE, toplabel="") {

  dataname <- deparse(substitute(data))
  outdir <- ms_directory(data, dataname)$outdir
  data <- ms_directory(data, dataname)$data

  if (length(complexID)) {
    if (class(complexID)=="character") { complexID <- as.numeric(complexID) }
    data <- subset(data, ComplexID %in% complexID)
    if (nrow(data)==0) { stop("Make sure the provided complexID is present in the input dataset") }
  }

  if (removeredundancy) {
    complexdata <- data %>% dplyr::group_by(ComplexID) %>%
      dplyr::summarise(subunitsIdentified=paste(id,collapse=";")) %>%
      dplyr::left_join(unique(data[ ,c("ComplexID","subunitsIdentifiedPerc")])) %>%
      dplyr::arrange(-subunitsIdentifiedPerc)

    i = 1
    while (i<nrow(complexdata)) {
      # l1 holds the members of query complex
      l1 <- stringr::str_split(as.character(complexdata[i,2]),";")[[1]]
      l <- c()
      for (k in (i+1):nrow(complexdata)) {
        # l2 holds the members of subject complex
        l2 <- stringr::str_split(as.character(complexdata[k,2]),";")[[1]]
        l <- c(l,sum(l1%in%l2)/min(length(l1),length(l2)))
      }
      if (length(which(l>=similaritythreshold))>0) {
        rem=-1*(i+which(l>=similaritythreshold))
        complexdata<-complexdata[rem, ]
      }
      #print(c(i,sum(l>=0.5),length(compdeffilt[,1])))
      i=i+1
    }
    data <- subset(data, ComplexID %in% complexdata$ComplexID)
    message(paste0("The number of multi-protein complexes after the removal of redundancy is: ",
            length(unique(data$ComplexID))))
  } else {
    message(paste0("The number of multi-protein complexes to analyze is: ",
                   length(unique(data$ComplexID))))
  }

  corr <- plyr::ddply(data, plyr::.(ComplexID), function(data) {
    if (sum(duplicated(data$gene))>0) {
      message("Found the following duplicated gene, need to be removed.")
      stop(message(data[duplicated(data$gene), c("ComplexID","id","gene")]))
    }
    name = unique(data[ ,2])
    genepos <- grep("^gene", names(data))
    data <- tidyr::gather(data[ ,c(genepos,7:(6+nread))], condition, reading, -gene)
    data1 <- tidyr::spread(data, gene, reading)
    correlation_table <- cor(data1[,-1], use="pairwise.complete.obs")
    #correlation_table1[lower.tri(correlation_table1, diag=T)] <- NA
    ut <- upper.tri(correlation_table)
    correlation_table <- data.frame(
      ComplexName = name,
      row = rownames(correlation_table)[row(correlation_table)[ut]],
      column = rownames(correlation_table)[col(correlation_table)[ut]],
      correlation = (correlation_table)[ut]
    )
  } )
  ms_filewrite(corr, paste0("correlation_in_", dataname, ".txt"),
               outdir=outdir, withdescription=F)

  corr_table <- corr %>% dplyr::group_by(ComplexID, ComplexName) %>%
    dplyr::summarise(mediancorr=median(correlation,na.rm=T), numofcorr=n(),
                  goodcorr=sum(correlation>goodcorrcutoff,na.rm=T),
                  vgoodcorr=sum(correlation>vgoodcorrcutoff,na.rm=T),
                  goodcorrperc=goodcorr/numofcorr, vgoodcorrperc=vgoodcorr/numofcorr) %>%
    dplyr::arrange(-mediancorr*log(numofcorr))
  ms_filewrite(corr_table, paste0("summarized_correlation_in_", dataname, ".txt"),
               outdir=outdir, withdescription=F)

  if (doplotting) {
    data$ComplexID <- factor(data$ComplexID, levels=corr_table$ComplexID)
    corr_plot <- plyr::dlply(data, plyr::.(ComplexID), function(data) {
      genepos <- grep("^gene", names(data))
      name = paste0("Complex ID: ", unique(data[,1]), "\n",
                    unique(data[,2]),"\n",
                    "Total subunit number: ", unique(data[,4]),"\n",
                    "Identified subunit number: ", nrow(data))
      data <- tidyr::gather(data[ ,c(genepos,7:(6+nread))], condition, reading, -gene)
      data1 <- tidyr::spread(data, gene, reading)
      q <- GGally::ggcorr(data1[,-1], palette="RdGy", size=3, label_round=2,
                          label=T, label_color="white", label_size=3, hjust=0.5)
      q <- q + ggtitle(name) +
        theme(
          text = element_text(size=8),
          plot.title = element_text(hjust=0.5, size=rel(1))
          #axis.text = element_text(angle=45, hjust=1, size=ref(0.7))
        )
    } )

    if (external) { external_graphs(T) }

    params <- list(nrow=layout[1], ncol=layout[2])
    n <- with(params, nrow*ncol)
    ## add one page if division is not complete
    pages <- length(corr_plot) %/% n + as.logical(length(corr_plot) %% n)
    groups <- split(seq_along(corr_plot), gl(pages, n, length(corr_plot)))

    pl <- lapply(names(groups), function(i) {
      gridExtra::grid.arrange(
        do.call(gridExtra::arrangeGrob,
                c(corr_plot[groups[[i]]], params, top=toplabel,left="",bottom="")))
    })

    class(pl) <- c("arrangelist", "ggplot", class(pl))
    ggplot2::ggsave(file=paste0(outdir,"/",format(Sys.time(), "%y%m%d_%H%M"),"_",
                                dataname,"_correlation_in_complex.pdf"), pl, height=11.69, width=8.27)

    if (external) { external_graphs(F) } # switch off the external graphs
  }
  if (length(attr(corr,"outdir"))==0 & length(outdir)>0) {
    attr(corr,"outdir") <- outdir
  }
  return(corr)
}
