#' imprints_imputation
#'
#' Function to impute the missing values in IMPRINTS-CETSA dataset.
#'
#'
#' @param data Normalized dataset but with missing values to be imputed and filled in
#' @param iteration an integer specifying how many rounds of imputation should be performed
#' @param percondition whether do the imputation on each condition separately, default
#' set to FALSE, to be implemented
#'
#' @import mice
#' @import VIM
#' @export
#' @return a dataframe
#' @examples \dontrun{
#'  MOLM <- imprints_imputation(MOLM)
#' }
#'
#'

imprints_imputation <- function(data, iteration=3, method="norm.predict",
                             percondition=FALSE, excludeNA=TRUE) {

  dataname <- deparse(substitute(data))
  outdir <- ms_directory(data, dataname)$outdir
  data <- ms_directory(data, dataname)$data

  if (length(grep("description", names(data)))) {
    proteininfo <- unique(data[ ,c("id","description")])
    data$description <- NULL
  }
  if (length(grep("countNum", names(data)))) {
    countinfo <- unique(data[ ,stringr::str_which(names(data), "^id$|^sumPSM|^countNum|^sumUniPeps")])
    data <- data[ ,-stringr::str_which(names(data), "^sumPSM|^countNum|^sumUniPeps")]
  }

  pdf(file = paste0(format(Sys.time(), "%y%m%d_%H%M_"),dataname,"missingness_pre.pdf"), width=14, height=6)
  mice::md.pattern(data[,-1], rotate.names=TRUE)
  dev.off()

  data1 <- data[,-1]
  names(data1) <- paste0("TEMP",names(data1))
  imp1 <- mice::mice(as.matrix(data1), m=iteration, maxit=5, method=method, printFlag=FALSE, seed=123)
  imp2 <- mice::complete(imp1, action="long", include=FALSE)
  names(imp2)[c(1:2)] <- c("imp","id")
  imp2$id <- data$id
  imp3 <- gather(imp2, condition, reading, -imp, -id)
  imp4 <- imp3 %>% dplyr::group_by(id,condition) %>% dplyr::summarise(reading=mean(reading,na.rm=T))
  imp4 <- tidyr::spread(imp4, condition, reading)
  names(imp4) <- gsub("TEMP","",names(imp4))

  pdf(file = paste0(format(Sys.time(), "%y%m%d_%H%M_"),dataname,"missingness_post.pdf"), width=14, height=6)
  mice::md.pattern(imp4[,-1], rotate.names=TRUE)
  dev.off()
  cat("Done...\n")

  if (excludeNA) {
    cat("To exclude any observations that still contain missing values after imputation...\n")
    imp4 <- na.exclude(imp4)
    pdf(file = paste0(format(Sys.time(), "%y%m%d_%H%M_"),dataname,"missingness_final.pdf"), width=14, height=6)
    mice::md.pattern(imp4[,-1], rotate.names=TRUE)
    dev.off()
    cat("Done...\n")
  }

  data1 <- merge(imp4, countinfo)
  data1 <- merge(proteininfo, data1)

  if (length(attr(data1,"outdir"))==0 & length(outdir)>0) {
    attr(data1,"outdir") <- outdir
  }
  ms_filewrite(data1, paste0(dataname,"_","IMPRINTS_imputed.txt"), outdir=outdir)

  return(data1)
}
