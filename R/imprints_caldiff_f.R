#' imprints_caldiff_f
#'
#' Function to calculate the pair-wise (per replicate and temperature) protein abundance differences,
#' This one returns the same results as imprints_caldiff(), but in a much faster manner
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
#'   MOLM <- imprints_caldiff_f(MOLM, reftreatment="DMSO")
#' }
#'
#'

imprints_caldiff_f <- function(data, reftreatment=NULL, withinrep=TRUE) {

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

  exp_design <- strsplit(colnames(data)[-1], "_")
  data <- as.data.frame(data)
  if (length(exp_design[[1]]) == 4) {
    set <- unique(unlist(lapply(exp_design, "[[", 1)))
    temperature <- unique(unlist(lapply(exp_design, "[[", 2)))
    replicate <- unique(unlist(lapply(exp_design, "[[", 3)))
    treatment <- colnames(data)[-grep(paste0("_", reftreatment, "$"), colnames(data))]
    if (withinrep) {
      for(t in temperature) {
        for(b in replicate) {
          ref <- grep(paste0(s, "_", t, "_", b, "_"), reftreatment_name, value = TRUE)
          treat <- grep(paste0(s, "_", t, "_", b, "_"), treatment, value = TRUE)
          data[ ,treat] <- data[ ,treat] - data[ ,ref]
        }
      }
    } else {
      for(s in set) {
        for(t in temperature) {
          ref <- grep(paste0(s, "_", t, "_"), reftreatment_name, value = TRUE)
          ref <- apply(data[ ,ref], 1, mean, na.rm = TRUE)
          treat <- grep(paste0(s, "_", t, "_"), treatment, value = TRUE)
          data[ ,treat] <- data[ ,treat] - ref
        }
      }
    }
    data[ ,reftreatment_name] <- 0
  } else if (length(exp_design[[1]]) == 3) {
    temperature <- unique(unlist(lapply(exp_design, "[[", 1)))
    replicate <- unique(unlist(lapply(exp_design, "[[", 2)))
    treatment <- colnames(data)[-grep(paste0("_", reftreatment, "$"), colnames(data))]
    if (withinrep) {
      for(t in temperature) {
        for(b in replicate) {
          ref <- grep(paste0(t, "_", b, "_"), reftreatment_name, value = TRUE)
          treat <- grep(paste0(t, "_", b, "_"), treatment, value = TRUE)
          data[ ,treat] <- data[ ,treat] - data[ ,ref]
        }
      }
    } else {
      for(t in temperature) {
        ref <- grep(paste0(t, "_"), reftreatment_name, value = TRUE)
        ref <- apply(data[ ,ref], 1, mean, na.rm = TRUE)
        treat <- grep(paste0(t, "_"), treatment, value = TRUE)
        data[ ,treat] <- data[ ,treat] - ref
      }
    }
    data[ ,reftreatment_name] <- 0
  } else {
    stop("make sure the namings of the columns of the dasaset are correct.")
  }

  data <- dplyr::as_tibble(data)
  data1 <- merge(data, countinfo)
  data1 <- merge(proteininfo, data1)

  if (length(attr(data1,"outdir"))==0 & length(outdir)>0) {
    attr(data1,"outdir") <- outdir
  }
  ms_filewrite(data1, paste0(dataname,"_","imprints_caldiff.txt"), outdir=outdir)
  return(data1)
}
