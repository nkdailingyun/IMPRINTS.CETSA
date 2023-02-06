#' imprints_average
#'
#' Function to calculate the averaged signals from the IMPRTINTS-CETSA result
#'
#' @param data dataset after calculating the relative protein abundance differences
#'
#' @import dplyr
#' @export
#' @return a dataframe
#' @examples \dontrun{
#'   MOLM <- imprints_average(MOLM)
#' }
#'
#'


imprints_average <- function(data) {

  dataname <- deparse(substitute(data))
  outdir <- ms_directory(data, dataname)$outdir
  data <- ms_directory(data, dataname)$data

  if (length(grep("description", names(data)))) {
    proteininfo <- unique(data[ ,c("id","description")])
    data$description <- NULL
  }
  if (length(grep("countNum", names(data)))) {
    countinfo <- unique(data[ ,c("id","sumUniPeps","sumPSMs","countNum")])
    data <- data[ ,!(names(data) %in% c("sumUniPeps","sumPSMs","countNum"))]
  }

  data <- gather(data, condition, reading, -id)
  if (length(unlist(strsplit(data$condition[1], "_")))==4) {
    data <- tidyr::separate(data, condition, into=c("set","temperature","replicate","treatment"), sep="_")
    data <- data %>% group_by(id,set,temperature,treatment) %>%
      summarise(reading.mean=mean(reading,na.rm=T))
    data <- tidyr::unite(data, condition, set, temperature, treatment, sep="_")
    data <- tidyr::spread(data, condition, reading.mean)
  } else if (length(unlist(strsplit(data$condition[1], "_")))==3) {
    data <- tidyr::separate(data, condition, into=c("temperature","replicate","treatment"), sep="_")
    data <- data %>% group_by(id,temperature,treatment) %>%
      summarise(reading.mean=mean(reading,na.rm=T))
    data <- tidyr::unite(data, condition, temperature, treatment, sep="_")
    data <- tidyr::spread(data, condition, reading.mean)
  } else {
    stop("make sure the namings of the columns of the dasaset are correct.")
  }

  data <- merge(data, countinfo)
  data <- merge(proteininfo, data)

  if (length(attr(data,"outdir"))==0 & length(outdir)>0) {
    attr(data,"outdir") <- outdir
  }
  return(data)
}
