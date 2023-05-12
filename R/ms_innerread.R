#' ms_innerread
#'
#' Internal function to parse and read in data from tab delimited files exported from Proteome Discoverer
#'
#' @param file file name to import
#' @param fchoose whether to choose file interactively
#' @param treatment names of treatments (temperature, dose, time) applied to samples
#' @param nread number of reading channels, should match the number of channels used
#' @param abdread whether to read in protein abundance data
#' @param PDversion which version of Proteome Discoverer the data is searched, possible values 20,21,22,24
#' @param fdrcontrol whether to check the protein FDR confidence level
#' @param refchannel names of reference channel used in Proteome Discoverer search, such as 126
#' @param channels names of the read-in channels
#' @param ... other arguments ignored (for compatibility with generic)
#'
#' @keywords internal
#'
#' @importFrom readr read_tsv locale spec date_names date_names_lang
#'
#' @return a dataframe
#' @examples \dontrun{
#'  ms_innerread("file.txt")
#' }
#'
#'
#'
ms_innerread <- function (file, fchoose, treatment, nread,
               fdrcontrol, refchannel, channels){
  if (length(treatment) != nread | length(channels) != nread) {
    stop("Make sure you specify the right number of channels!")
  }
  if (length(file) == 0) {
    stop("No valid file(s) loaded!")
  }
  if (fchoose) {
    data <- readr::read_tsv(file.choose(), show_col_types = FALSE)
  }
  else {
    data <- readr::read_tsv(file = file, show_col_types = FALSE)
  }

  conditions <- grep(pattern = "^Abundance: ",
                     colnames(data), value = TRUE)
  conditions <- gsub(pattern = ".*, ", "", conditions)
  conditions <- unique(conditions)
  if (length(conditions) != 1) {
    stop("Make sure the condition was correctedly specified. Specify one unique condition for each input file.")
  }

  names(data) <- gsub(pattern = "_", "", names(data))
  names(data) <- gsub(pattern = "\\(", "", names(data))
  names(data) <- gsub(pattern = "\\)", "", names(data))

  if (fdrcontrol) {
    pattern <- grep("Protein FDR Confidence", names(data),
                    value = FALSE)
    if (length(pattern) > 0) {
      names(data)[pattern] <- "FDR"
      pattern <- grep("Low", data$FDR)
      if (length(pattern) > 0) {
        message(paste0("There are ", length(pattern),
                       " low confidence proteins in the origial data."))
        message("These low confidence proteins were removed from downstream analysis!")
        data_lc <- data[pattern, ]
        data <- data[-pattern, ]
      }
      else {
        message("No low confidence proteins found in the origial data.")
      }
    }
  }
  colnames <- names(data)

  pattern <- grep("Standard Error", colnames)
  if (length(pattern) > 0) {
    data <- data[, -pattern]
    colnames <- names(data)
  }
  pattern <- grep("Variability", colnames)
  if (length(pattern) > 0) {
    data <- data[, -pattern]
    colnames <- names(data)
  }
  pattern <- grep("Grouped CV", colnames)
  if (length(pattern) > 0) {
    data <- data[, -pattern]
    colnames <- names(data)
  }

  id <- grep(pattern = "^Accession$", colnames)
  description <- grep(pattern = "^Description$", colnames)
  sumUnipeps <-  grep(pattern = "^# Unique Peptides$", colnames)
  sumPSMs <-  grep(pattern = "^# PSMs$", colnames)
  countNum <- grep(pattern = paste0("^Abundances Count [A-z0-9,. -]+",
                        refchannel, "[A-z0-9,. -]+$"), colnames)
  raw_abundance <- grep(pattern = "^Abundance: ", colnames)

  data <- data[,c(id, description, raw_abundance, sumUnipeps, sumPSMs, countNum)]

  if(length(countNum) == 0){
    message("There is no protein abundance count information in the original data.")
    message("Set a default number for the protein abundance count to 2.")
    data$countNum <- 2
  }

  data$condition <- conditions
  colnames(data) <- c("id", "description",
                      treatment, "sumUnipeps",
                      "sumPSMs", "countNum", "condition")
  data <- data[,c(1:2, ncol(data), 3:(ncol(data) - 1))]


  return(data)
}
