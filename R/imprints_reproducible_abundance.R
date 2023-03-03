#' imprints_reproducible_abundance
#'
#' Function to retrieve a reproducible subset of the IMPRINTS-CETSA dataset,
#' regarding to both expression level changes and thermal shifts.
#' This is different from imprints_reproducible() function,
#' which is performed on the logratio data but not normalized dataset with protein abundance.
#'
#' @param data dataset after imprints_normalization() function
#' @param set a single character to indicate the sample name
#' @param treatment a single character to indicate the sample name
#' @param cvthreshold the CV threshold value for subsetting reproducible
#' measurements, default value is 0.1
#' @param returnlist whether to return as list but not as dataframe,
#' default set to FALSE
#'
#' @import ggpubr
#' @import ggrepel
#' @export
#' @return a list of dataframe, for each condition
#' @examples \dontrun{
#'   MOLM2 <- imprints_reproducible_abundance(MOLM1)
#' }
#'

imprints_reproducible_abundance <- function(data, set=NULL, treatment=NULL,
                               cvthreshold=0.1, returnlist=FALSE) {

  dataname <- deparse(substitute(data))
  outdir <- ms_directory(data, dataname)$outdir
  data <- ms_directory(data, dataname)$data
  data_copy <- data

  if (length(grep("description", names(data)))) {
    proteininfo <- unique(data[ ,c("id","description")])
  }
  if (length(grep("countNum", names(data)))) {
    countinfo <- unique(data[ ,stringr::str_which(names(data), "^id$|^sumPSM|^countNum|^sumUniPeps")])
    data <- data[ ,-stringr::str_which(names(data), "^sumPSM|^countNum|^sumUniPeps")]
  }

  if (length(set)==1) { data <- data[ ,c(1,2,grep(set,names(data)))] }
  else if (length(set)>1) { stop("Please specify only one set name") }

  datal <- tidyr::gather(data[ ,-2], treatment, reading, -id)
  if (length(set)>0) {
    if(!length(unlist(strsplit(datal$treatment[1], "_")))==4) {
      stop("Please make sure there is set information incorporated in the column head")
    }
    datal1 <- tidyr::separate(datal, treatment, into=c("set","temperature","replicate","treatment"), sep="_")
    datal2 <- list()
    for (i in unique(datal1$treatment)) {
      datal_temp <- subset(datal1, treatment==i)
      message(paste0("in condition: ", i))
      datal_score <- datal_temp %>% dplyr::group_by(id, set, treatment, temperature) %>%
        dplyr::summarise(mreading=mean(reading,na.rm=T), sdreading=sd(reading, na.rm=T)) %>%
        dplyr::group_by(id) %>% dplyr::mutate(cvreading=sqrt(exp(sdreading^2)-1)) %>%
        dplyr::group_by(id) %>% dplyr::mutate(mcvreading=mean(cvreading,na.rm=T))
      write.csv(datal_score, paste0(outdir,"/",format(Sys.time(),"%y%m%d_%H%M_"),i,"_abundance_CV.csv"), row.names=F)
      datal_score <- datal_score %>% dplyr::group_by(id) %>% summarise(mcv=mean(mcvreading,na.rm=T))
      q <- ggpubr::ggdensity(datal_score, x="mcv", fill="lightgray", add="median", rug=F, xlab="CV", title=as.character(i))
      ggplot2::ggsave(file=paste0(outdir,"/",format(Sys.time(), "%y%m%d_%H%M"), "_", i, "_abundance_CV.pdf"), q, width=4, height=4)
      reproducible1 <- datal_score %>% filter(mcv<=cvthreshold)
      message(paste0(nrow(reproducible1), " proteins passed the cv cutoff."))
      message(paste0("The median + 2*MAD cutoff of CV is ", round(median(datal_score$mcv,na.rm=T)+2*mad(datal_score$mcv,na.rm=T),3)))
      datal_temp <- subset(datal_temp, id %in% reproducible1$id)
      datal_temp <- tidyr::unite(datal_temp, condition, temperature, replicate, treatment, sep="_")
      datal_temp <- tidyr::spread(datal_temp, condition, reading)
      datal_temp <- merge(datal_temp, countinfo)
      datal_temp <- merge(proteininfo, datal_temp)
      if (length(attr(datal_temp,"outdir"))==0 & length(outdir)>0) {
        attr(datal_temp,"outdir") <- outdir
      }
      datal2[[i]] <- datal_temp
    }
  } else {
    if(!length(unlist(strsplit(datal$treatment[1], "_")))==3) {
      stop("Please make sure there is treatment information incorporated in the column head")
    }
    datal1 <- tidyr::separate(datal, treatment, into=c("temperature","replicate","treatment"), sep="_")
    datal2 <- list()
    for (i in unique(datal1$treatment)) {
      datal_temp <- subset(datal1, treatment==i)
      message(paste0("in condition: ", i))
      datal_score <- datal_temp %>% dplyr::group_by(id, treatment, temperature) %>%
        dplyr::summarise(mreading=mean(reading,na.rm=T), sdreading=sd(reading, na.rm=T)) %>%
        dplyr::group_by(id) %>% dplyr::mutate(cvreading=sqrt(exp(sdreading^2)-1)) %>%
        dplyr::group_by(id) %>% mutate(mcvreading=mean(cvreading,na.rm=T))
      write.csv(datal_score, paste0(outdir,"/",format(Sys.time(),"%y%m%d_%H%M_"),i,"_abundance_CV.csv"), row.names=F)
      datal_score <- datal_score %>% group_by(id) %>% summarise(mcv=mean(mcvreading,na.rm=T))
      q <- ggpubr::ggdensity(datal_score, x="mcv", fill="lightgray", add="median", rug=F, xlab="CV", title=as.character(i))
      ggplot2::ggsave(file=paste0(outdir,"/",format(Sys.time(), "%y%m%d_%H%M"), "_", i, "_abundance_CV.pdf"), q, width=4, height=4)
      reproducible1 <- datal_score %>% dplyr::filter(mcv<=cvthreshold)
      message(paste0(nrow(reproducible1), " proteins passed the cv cutoff."))
      message(paste0("The median + 2*MAD cutoff of CV is ", round(median(datal_score$mcv,na.rm=T)+2*mad(datal_score$mcv,na.rm=T),3)))
      datal_temp <- subset(datal_temp, id %in% reproducible1$id)
      datal_temp <- tidyr::unite(datal_temp, condition, temperature, replicate, treatment, sep="_")
      datal_temp <- tidyr::spread(datal_temp, condition, reading)
      datal_temp <- merge(datal_temp, countinfo)
      datal_temp <- merge(proteininfo, datal_temp)
      if (length(attr(datal_temp,"outdir"))==0 & length(outdir)>0) {
        attr(datal_temp,"outdir") <- outdir
      }
      datal2[[i]] <- datal_temp
    }
  }
  if (returnlist) {
    return(datal2)
  } else {
    repid <- NULL
    for (i in names(datal2)) {
      repid <- unique(c(repid, datal2[[i]]$id))
    }
    data <- subset(data_copy, id %in% repid)
    if (length(attr(data,"outdir"))==0 & length(outdir)>0) {
      attr(data,"outdir") <- outdir
    }
    return(data)
  }
}
