#' ms_innerread
#'
#' Internal function to parse and read in data from tab delimited files exported from Proteome Discoverer
#'
#' @param file file name to import
#' @param fchoose whether to choose file interactively
#' @param treatment names of treatments (temperature, dose, time) applied to samples
#' @param nread number of reading channels, should match the number of channels used
#' @param fdrcontrol whether to check the protein FDR confidence level
#' @param refchannel names of reference channel used in Proteome Discoverer search, such as 126
#' @param channels names of the read-in channels
#'
#' @keywords internal
#'
#' @importFrom readr read_tsv locale spec date_names date_names_lang
#'
#' @return a dataframe
#'
#'
ms_innerread <- function(file, fchoose, treatment, nread,
                         fdrcontrol, refchannel, channels) {

  # Treatment and channel check
  if (length(treatment) != nread | length(channels) != nread) {
    stop("Make sure you specify the right number of channels!")
  }
  # File Check
  if (length(file) == 0) {
    stop("No valid file(s) loaded!")
  }
  # File Loading
  if (fchoose) {
    data <- readr::read_tsv(file.choose(),show_col_types=FALSE)
  } else {
    data <- readr::read_tsv(file=file, show_col_types=FALSE)
  }

  # Clean up column names
  names(data) <- gsub(pattern="_", "", names(data))
  names(data) <- gsub(pattern="\\(", "", names(data))
  names(data) <- gsub(pattern="\\)", "", names(data))

  # create condition vector
  conditions <- grep(pattern = "^Abundance: |^Abundances: ", colnames(data), value = TRUE)
  if (length(conditions)==0) {
    conditions <- grep(pattern = "^Abundances Grouped: ", colnames(data), value = TRUE)
  }
  conditions <- gsub(pattern = ".*, ", "", conditions)
  conditions <- unique(conditions)
  if (length(conditions) != 1) {
    stop("Make sure the condition was correctedly specified. Specify one unique condition for each input file.")
  }
  # print(conditions)

  # Remove low confidence proteins
  if (fdrcontrol) {
    pattern <- grep("Protein FDR Confidence", names(data), value=FALSE)
    if (length(pattern) > 0) {
      names(data)[pattern] <- "FDR"
      pattern <- grep("Low", data$FDR)
      if (length(pattern) > 0) {
        message(paste0("There are ",length(pattern)," low confidence proteins in the origial data."))
        message("These low confidence proteins are removed from downstream analysis!")
        data_lc <- data[pattern, ]
        data <- data[-pattern, ]
      } else {
        message("No low confidence proteins found in the origial data.")
      }
    }
  }

  colnames <- names(data)

  # remove standard error and variability
  pattern <- grep("Standard Error", names(data), value=FALSE)
  if (length(pattern) > 0) {
    data <- data[ ,-pattern]
    colnames <- names(data)
  }
  pattern <- grep("Variability", names(data), value=FALSE)
  if (length(pattern) > 0) {
    data <- data[ ,-pattern]
    colnames <- names(data)
  }
  pattern <- grep("Grouped CV", names(data), value=FALSE)
  if (length(pattern) > 0) {
    data <- data[ ,-pattern]
    colnames <- names(data)
  }

  id <- grep(pattern = "^Accession$", colnames)
  description <- grep(pattern = "^Description$", colnames)
  sumUnipeps <-  grep(pattern = "^# Unique Peptides$", colnames)
  sumPSMs <-  grep(pattern = "^# PSMs$", colnames)
  countNum <- grep(pattern = paste0("^Abundances Count: [A-z0-9,. -]+",": ",
                                    refchannel, "[A-z0-9,. -]+$"), colnames)
  raw_abundance <- grep(pattern = "^Abundance: |^Abundances: ", colnames)
  if (length(raw_abundance)==0) {
    raw_abundance <- grep(pattern = "^Abundances Grouped: ", colnames)
  }

  data <- data[ ,c(id, description, raw_abundance, sumUnipeps, sumPSMs, countNum)]

  if (length(countNum) == 0) {
    message("There is no protein abundance count information in the original data.")
    message("Set a default number for the protein abundance count to 2.")
    data$countNum <- 2
  }
  colnames(data)[c((ncol(data)-3):(ncol(data)-1))] <- c("sumUniPeps", "sumPSMs", "countNum")

  data$condition <- conditions

  # be sure that channels are in same order as treatment vector
  ord_channel <- sapply(channels,
                        function(y) grep(paste0(" ", y, ","), colnames(data)))
  data <- data[,c(1:2, ord_channel, (ncol(data)-3):ncol(data))]

  # rename columns
  colnames(data) <- c("id", "description", treatment, "sumUniPeps",
                      "sumPSMs", "countNum", "condition")
  data <- data[,c(1:2, ncol(data), 3:(ncol(data)-1))]

  return(data)
}
