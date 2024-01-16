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
#' @param software Software used to obtain Protein Group files. Either PD, MaxQuant or pFind.
#'   Default is PD
#'
#' @keywords internal
#'
#' @importFrom readr read_tsv locale spec date_names date_names_lang
#'
#' @return a dataframe
#'
#'
ms_innerread <- function(file, fchoose, treatment, nread,
                         fdrcontrol, refchannel, channels, software) {

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
  if(software == "PD"){
    conditions <- grep(pattern = "^Abundance: |^Abundances: ", colnames(data), value = TRUE)
    if (length(conditions)==0) {
      conditions <- grep(pattern = "^Abundances Grouped: ", colnames(data), value = TRUE)
    }
    conditions <- gsub(pattern = ".*, ", "", conditions)
    conditions <- unique(conditions)
    if (length(conditions) != 1) {
      stop("Make sure the condition was correctedly specified. Specify one unique condition for each input file.")
    }
  }
  else if(software == "MaxQuant"){
    conditions <- grep(pattern = "^Reporter intensity corrected \\d{1,}", colnames(data), value = TRUE)
    if (length(conditions)==0) {
      conditions <- grep(pattern = "^Reporter intensity \\d{1,}", colnames(data), value = TRUE)
    }
    conditions <- gsub(pattern = ".* ", "", conditions)
    conditions <- unique(conditions)
    if (length(conditions) != 1) {
      stop("Make sure the condition was correctedly specified. Specify one unique condition for each input file.")
    }
  }
  # print(conditions)

  # Remove low confidence proteins
  if (fdrcontrol) {
    if(software == "PD"){
      pattern <- grep("Protein FDR Confidence", names(data), value=FALSE)
    }
    else if(software == "MaxQuant"){
      pattern <- grep("^Q-value$", names(data), value=FALSE)
    }

    if (length(pattern) == 1) {
      if(software == "PD"){
        pattern <- grep("Low", data[[pattern]])
      }
      else if(software == "MaxQuant"){
        pattern <- which(data[[pattern]] <= 0)
      }

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

  if(software == "PD"){
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
  }
  else if(software == "MaxQuant"){
    # removing proteins marked as reversed if any
    pattern <- grep("^(R|r)everse$", names(data), value=FALSE)
    if (length(pattern) > 0) {
      reversed <- which(data[[pattern]] == "+")
      if(length(reversed)){
        data <- data[-reversed,]
      }
    }

    id <- grep(pattern = "^Majority protein IDs$", colnames)
    description <- grep(pattern = "^Fasta (H|h)eaders$", colnames)
    sumUnipeps <-  grep(pattern = "^Unique (P|p)eptides$", colnames)
    sumPSMs <-  grep(pattern = "^MS/MS (c|C)ount$", colnames)
    countNum <- grep(pattern = paste0("^Reporter intensity count ",
                                      which(channels == refchannel), " "), colnames)

    raw_abundance <- grep(pattern = "^Reporter intensity corrected \\d{1,}", colnames(data))
    if (length(raw_abundance)==0) {
      raw_abundance <- grep(pattern = "^Reporter intensity \\d{1,}", colnames(data))
    }
  }

  data <- data[ ,c(id, description, raw_abundance, sumUnipeps, sumPSMs, countNum)]

  if (length(countNum) == 0) {
    message("There is no protein abundance count information in the original data.")
    message("Set a default number for the protein abundance count to 2.")
    data$countNum <- 2
  }

  data$condition <- conditions
  colnames(data)[c((ncol(data)-3):(ncol(data)-1))] <- c("sumUniPeps", "sumPSMs", "countNum")

  if(software == "PD"){
    # be sure that channels are in same order as treatment vector
    ord_channel <- sapply(channels,
                          function(y) grep(paste0(" ", y, ","), colnames(data)))
    data <- data[,c(1:2, ord_channel, (ncol(data)-3):ncol(data))]
  }

  # rename columns
  colnames(data) <- c("id", "description", treatment, "sumUniPeps",
                      "sumPSMs", "countNum", "condition")
  data <- data[,c(1:2, ncol(data), 3:(ncol(data)-1))]

  if(software == "MaxQuant"){
    # reformat description
    data$description <- unlist(unname(sapply(data$description,
                                             function(x){
                                               x <- strsplit(x, ";")[[1]]
                                               x <- unlist(lapply(strsplit(x, "\\|"),
                                                                  function(y){
                                                                    if(length(y) < 3){
                                                                      y <- NULL
                                                                    }
                                                                    else{
                                                                      y <- y[[3]]
                                                                    };
                                                                    y
                                                                  })
                                               )
                                               x <- sub(".*? ", "", x)
                                               x <- paste(x, collapse = ";");
                                               x
                                             }, simplify = FALSE)
                                      )
                               )
  }

  return(data)
}
