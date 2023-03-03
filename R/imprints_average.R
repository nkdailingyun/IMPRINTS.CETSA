#' imprints_average
#'
#' Function to calculate the averaged signals from the IMPRINTS-CETSA result
#'
#' @param data dataset after calculating the relative protein abundance differences
#' @param savefile a logical to tell if you want save the results or not
#'
#' @importFrom tidyr expand gather separate unite
#' @importFrom dplyr group_by summarise ungroup
#' @export
#' @return a dataframe
#' @examples \dontrun{
#'   MOLM <- imprints_average(MOLM)
#' }
#'
#'

imprints_average <- function(data, savefile=TRUE) {

  if (savefile) {
    filename <- paste0(deparse(substitute(data)), "_average", ".txt")
  }
  dataname <- deparse(substitute(data))
  outdir <- ms_directory(data, dataname)$outdir
  data <- ms_directory(data, dataname)$data

  if (length(grep("description", names(data)))) {
    proteininfo <- unique(data[ ,c("id","description")])
    data$description <- NULL
  }
  # if (length(grep("countNum", names(data)))) {
  #   countinfo <- unique(data[ ,c("id","sumUniPeps","sumPSMs","countNum")])
  #   data <- data[ ,!(names(data) %in% c("sumUniPeps","sumPSMs","countNum"))]
  # }
  if (length(grep("countNum", names(data)))) {
    countinfo <- unique(data[ ,stringr::str_which(names(data), "^id$|^sumPSM|^countNum|^sumUniPeps")])
    data <- data[ ,-stringr::str_which(names(data), "^sumPSM|^countNum|^sumUniPeps")]
    #allows to work with joined table
  }

  data <- tidyr::gather(data, condition, reading, -id)
  if (length(unlist(strsplit(data$condition[1], "_")))==4) {
    data <- tidyr::separate(data, condition, into=c("set","temperature","replicate","treatment"), sep="_")
    data <- data %>% dplyr::group_by(id,set,temperature,treatment) %>%
      dplyr::summarise(reading.mean=mean(reading,na.rm=T))
    data <- tidyr::unite(data, condition, set, temperature, treatment, sep="_")
    data <- tidyr::spread(data, condition, reading.mean)
  } else if (length(unlist(strsplit(data$condition[1], "_")))==3) {
    data <- tidyr::separate(data, condition, into=c("temperature","replicate","treatment"), sep="_")
    data <- data %>% dplyr::group_by(id,temperature,treatment) %>%
      dplyr::summarise(reading.mean=mean(reading,na.rm=T))
    data <- tidyr::unite(data, condition, temperature, treatment, sep="_")
    data <- tidyr::spread(data, condition, reading.mean)
  } else {
    stop("make sure the namings of the columns of the dasaset are correct.")
  }

  data <- merge(data, countinfo)
  data <- merge(proteininfo, data)

  if (savefile) {
    ms_filewrite(data, filename)
  }

  if (length(attr(data,"outdir"))==0 & length(outdir)>0) {
    attr(data,"outdir") <- outdir
  }
  return(data)
}
