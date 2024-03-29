#' imprints_rearrange
#'
#' Function to rearrange IMPRINTS-CETSA data into a wide one-unique-protein-per-row format
#'
#' @param data dataset to transform into a wide format
#' @param nread number of reading channels, default value 9 in typical IMPRINTS-CETSA scheme
#' @param repthreshold the minimal percentage threshold of protein being sampled
#' from multiple runs, default value is 0.75
#' @param withabdreading whether the kept proteins should have readings at abundance
#' reference channel, usually 37C or the lowest heating temperature (basetemp)
#' @param basetemp the character indicating the baseline temperature for
#' expression levels, default value is 37C
#' @param averagecount whether to take the median of the abundance count numbers across
#' the measured temperature range and then use this value for filtering, default set to TRUE,
#' otherwise, filter the proteins according to the associated count numbers at each temperature
#' @param countthreshold the minimal threshold number of associated abundance
#' count of proteins, default value is 2
#' @param extraid a vector of UniprotID to be included in the subset dataset
#'
#' @importFrom tidyr expand gather separate unite
#' @importFrom dplyr filter group_by inner_join left_join mutate n rowwise summarise top_n ungroup
#' @export
#' @return a dataframe
#' @examples \dontrun{
#'  MOLM <- imprints_rearrange(MOLM, nread=9, countthreshold=2)
#' }
#'
#'


imprints_rearrange <- function(data, nread=9, repthreshold=0.75, withabdreading=TRUE,
                               basetemp="37C", averagecount=TRUE, countthreshold=2,
                               extraid=NULL) {

  dataname <- deparse(substitute(data))
  outdir <- ms_directory(data, dataname)$outdir
  data <- ms_directory(data, dataname)$data
  data_copy <- data

  if (!averagecount & countthreshold>0) { # filter on each temperature
    data <- subset(data, countNum >= countthreshold)

    message(paste0(length(unique(data$id)), " proteins pass the cutoff of PSM count number ",
                   countthreshold, "."))
  }

  if (repthreshold>0) {
    counttable <- dplyr::count(data, id) %>%
      dplyr::mutate(freq = n/max(n)) %>%
      dplyr::filter(freq >= repthreshold)
    data <- subset(data, id %in% counttable$id)

    message(paste0(nrow(counttable), " proteins pass the measurement replicates cutoff ",
                   repthreshold*100, "%."))
  }

  if (withabdreading) {
    refcond <- unique(data$condition)[grep(basetemp, unique(data$condition))]

    if (length(refcond)==1) {
      refdata <- subset(data, condition == refcond)
    }
    else if (length(refcond)>1) {
      refdata <- subset(data, condition %in% refcond)
      refdata <- dplyr::count(refdata,id) %>%
        dplyr::filter(n==length(refcond))
    }

    message(paste0(length(unique(refdata$id)), " proteins are measured at ",
                   basetemp," and they are kept."))
    data <- subset(data, id %in% refdata$id)
  }

  if (length(extraid)) {
    message(paste0("To include extra ", length(extraid), " proteins as specified."))
    extraid <- setdiff(extraid, unique(data$id))
    data <- rbind(data, data_copy[which(data_copy$id %in% extraid), ])
  }

  d1 <- tidyr::gather(data[ ,c(1,3,4:(3+nread))], treatment, reading, -id, -condition)
  d1 <- tidyr::unite(d1, combinedcol, condition, treatment, sep="_")
  d1 <- tidyr::spread(d1, combinedcol, reading)
  data1 <- unique(data[ ,c(1,2)]) %>% dplyr::inner_join(d1, by="id")
  #%>% rowwise() %>% mutate(gene=getGeneName(description))

  # extract median PSM, peptide and count numbers
  counttable <- NULL
  peppos <- grep("^sumUniPeps", names(data), value=F)
  psmpos <- grep("^sumPSMs", names(data), value=F)
  countpos <- grep("^countNum", names(data), value=F)
  pos <- c(peppos, psmpos, countpos)
  if (length(pos) > 1) {
    counttable <- dplyr::group_by(data[ ,c(1,pos)], id) %>%
      dplyr::summarise(sumUniPeps = median(sumUniPeps),
                       sumPSMs = median(sumPSMs),
                       countNum = median(countNum))

    if (averagecount & countthreshold > 0) { # filter on median across temperatures
      counttable <- subset(counttable, countNum>=countthreshold)
      message(paste0(nrow(counttable), " proteins pass the cutoff of PSM count number ",
                     countthreshold, "."))
    }
    data1 <- merge(data1, counttable)
  }

  if (length(attr(data1,"outdir"))==0 & length(outdir)>0) {
    attr(data1,"outdir") <- outdir
  }
  ms_filewrite(data1, paste0(dataname, "_data_pre_normalization.txt"),
                          outdir=outdir)
  return(data1)
}
