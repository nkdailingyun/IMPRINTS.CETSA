#' imprints_caldiff
#'
#' Function to calculate the pair-wise (per replicate and temperature) protein abundance differences
#'
#' @param data Normalized dataset to calculate the relative protein abundance differences
#' @param reftreatment the name of treatment used as the reference control, such as "DMSO" or "Ctrl"
#' @param withinrep whether the calculation of the relative protein abundance difference should
#' still within the same biorep, default set to TRUE, when the bioreps are balanced
#'
#' @importFrom tidyr expand gather separate unite
#' @export
#' @return a dataframe
#' @examples \dontrun{
#'   MOLM <- imprints_caldiff(MOLM, reftreatment="DMSO")
#' }
#'
#'

imprints_caldiff <- function(data, reftreatment=NULL, withinrep=TRUE) {

  options(digits = 5)
  dataname <- deparse(substitute(data))
  outdir <- ms_directory(data, dataname)$outdir
  data <- ms_directory(data, dataname)$data

  if (length(reftreatment)==0) {
    stop("Need to specify the reference condition")
  }
  reftreatment_name <- grep(paste0("_", reftreatment, "$"), colnames(data), value = TRUE)
  if (!length(reftreatment_name)) {
    stop("Need to specify a right treatment condition as reference")
  }
  if (length(grep("description", names(data)))) {
    proteininfo <- unique(data[ ,c("id","description")])
    data$description <- NULL
  }
  if (length(grep("countNum", names(data)))) {
    countinfo <- unique(data[ ,stringr::str_which(names(data), "^id$|^sumPSM|^countNum|^sumUniPeps")])
    data <- data[ ,-stringr::str_which(names(data), "^sumPSM|^countNum|^sumUniPeps")]
  }

  data <- tidyr::gather(data, condition, reading, -id)
  if (length(unlist(strsplit(data$condition[1], "_")))==4) {
    data <- tidyr::separate(data, condition, into=c("set","temperature","replicate","treatment"), sep="_")
    if (sum(grepl(reftreatment, sort(unique(data$treatment))))!=1) {
      stop("Need to specify a right treatment condition as reference")
    }
    treatmentlevel <- c(reftreatment,setdiff(unique(data$treatment),reftreatment))
    data$treatment <- factor(data$treatment, levels=treatmentlevel)
    if (withinrep) {
      data1 <- plyr::ddply(data, c("set","temperature","replicate","id"), function(data) {
        data<-data[order(data$treatment), ]
        base=data$reading[1]
        data<-mutate(data, reading=reading-base)
      })
    } else {
      data1 <- plyr::ddply(data1, c("set","temperature","id"), function(data) {
        data<-data[order(data$treatment), ]
        base=mean(subset(data, treatment==reftreatment)$reading,na.rm=T)
        data<-mutate(data, reading=ifelse(treatment==reftreatment, 0.0, reading-base))
      })
    }
    data1 <- tidyr::unite(data1, condition, set, temperature, replicate, treatment, sep="_")
  } else if (length(unlist(strsplit(data$condition[1], "_")))==3) {
    data <- tidyr::separate(data, condition, into=c("temperature","replicate","treatment"), sep="_")
    if(sum(grepl(reftreatment, sort(unique(data$treatment))))!=1) {
      stop("Need to specify a right treatment condition as reference")
    }
    treatmentlevel <- c(reftreatment,setdiff(unique(data$treatment),reftreatment))
    data$treatment <- factor(data$treatment, levels=treatmentlevel)
    if (withinrep) {
      data1 <- plyr::ddply(data, c("temperature","replicate","id"), function(data) {
        data<-data[order(data$treatment), ]
        base=data$reading[1]
        data<-mutate(data, reading=reading-base)
      })
    } else {
      data1 <- plyr::ddply(data, c("temperature","id"), function(data) {
        data<-data[order(data$treatment), ]
        base=mean(subset(data, treatment==reftreatment)$reading,na.rm=T)
        data<-mutate(data, reading=ifelse(treatment==reftreatment, 0.0, reading-base))
      })
    }
    data1 <- tidyr::unite(data1, condition, temperature, replicate, treatment, sep="_")
  } else {
    stop("make sure the namings of the columns of the dasaset are correct.")
  }

  data1 <- tidyr::spread(data1, condition, reading)
  data1 <- merge(data1, countinfo)
  data1 <- merge(proteininfo, data1)

  if (length(attr(data1,"outdir"))==0 & length(outdir)>0) {
    attr(data1,"outdir") <- outdir
  }
  ms_filewrite(data1, paste0(dataname,"_","imprints_caldiff.txt"), outdir=outdir)
  return(data1)
}
