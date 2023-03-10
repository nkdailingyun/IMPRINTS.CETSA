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
ms_innerread <- function(file, fchoose, treatment, nread, abdread,
                         PDversion, fdrcontrol, refchannel, channels) {

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
  names(data) <- gsub(pattern=":", "", names(data))
  names(data) <- gsub(pattern="\\(", "", names(data))
  names(data) <- gsub(pattern="\\)", "", names(data))

  # Remove low confidence proteins
  if (fdrcontrol) {
    pattern <- grep("Protein FDR Confidence", names(data), value=FALSE)
    if (length(pattern) > 0) {
      names(data)[pattern] <- "FDR"
      pattern <- grep("Low", data$FDR)
      if (length(pattern) > 0) {
        message(paste0("There are ",length(pattern)," low confidence proteins in the origial data."))
        message("These low confidence proteins were removed from downstream analysis!")
        data_lc <- data[pattern, ]
        data <- data[-pattern, ]
      } else {
        message("No low confidence proteins found in the origial data.")
      }
    }
  }

  nrowdata <- nrow(data)
  colnames <- names(data)
  collength <- length(names(data))

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

  # create condition vector
  #conditionchannel <- setdiff(channels, refchannel)[1]
  # create condition vector
  conditions <- grep(pattern=paste0("^[A-z0-9,. -]+ / ", refchannel, "[A-z0-9,. -]+$"), colnames, value=TRUE)
  conditions <- gsub(pattern=paste0("^[A-z0-9,. -]+ / ", refchannel,", "), "", conditions)
  conditions <- unique(conditions)
  # print(conditions)

  # Move Accession
  tmppos <- grep(pattern="^Accession$", names(data), value=FALSE)
  collength <- length(names(data))
  data <- data[ ,c(tmppos,1:(tmppos-1),(tmppos+1):collength)]

  # Move Description
  tmppos <- grep(pattern="^Description$", names(data), value=FALSE)
  collength <- length(names(data))
  data <- data[ ,c(1,tmppos,2:(tmppos-1),(tmppos+1):collength)]

  # Condition row
  if (length(conditions) == 1) {
    data$condition <- rep(conditions, nrowdata)
    collength <- length(names(data))
    data <- data[ ,c(1:2,collength,3:(collength-1))]
  } else {
    stop("Make sure the condition was correctedly specified. Specify one unique condition for each input file.")
  }

  # Move reading values and rename them to correct position
  j = 1
  for (i in seq_len(nread)) {
    tmppos <- grep(pattern=paste0(channels[i],", ",conditions," / ",refchannel,", ",conditions), names(data), value=FALSE)
    if (length(tmppos)) {
      collength <- length(names(data))
      data <- data[ ,c(1:(j+2),tmppos,(j+3):(tmppos-1),(tmppos+1):collength)]
      j = j + 1
    } else if (channels[i] == refchannel) {
      # Set up start reference channel
      data$Reference <- rep(1.0, nrowdata)
      collength <- length(names(data))
      data <- data[ ,c(1:(j+2),collength,(j+3):(collength-1))]
      j = j + 1
    }
  }

  # to read in protein abundance raw data
  if (PDversion>=21 & abdread) { # PD2.0 doesn't contain abundance data
    for (i in seq_len(nread)) {
      #get column for each channel and move it to correct position:
      tmppos1 <- grep(pattern=paste0("^Abundances Grouped ",channels[i],"[A-z0-9,. -]+$"), names(data), value=FALSE)
      tmppos2 <- grep(pattern=paste0("^Abundance F[0-9]+ ",channels[i],"[A-z0-9,. -]+$"), names(data), value=FALSE)
      tmppos <- unique(c(tmppos1,tmppos2))
      collength <- length(names(data))
      data <- data[ ,c(1:(nread+2+i),tmppos,(nread+3+i):(tmppos-1),(tmppos+1):collength)]
    }
  }

  # Unique Peptides & PSMs
  tmppos <- grep(pattern="^# Unique Peptides$", names(data), value=FALSE)
  collength <- length(names(data))
  if (PDversion>=21 & abdread) {
    data <- data[,c(1:(2*nread+3),tmppos,(2*nread+4):collength)]
  } else {
    data <- data[,c(1:(nread+3),tmppos,(nread+4):collength)]
  }

  tmppos <- grep(pattern="^# PSMs$", names(data), value=FALSE)
  collength <- length(names(data))
  if (PDversion>=21 & abdread) {
    data <- data[ ,c(1:(2*nread+4),tmppos,(2*nread+5):collength)]
  } else {
    data <- data[ ,c(1:(nread+4),tmppos,(nread+5):collength)]
  }

  # Abundance counts
  tmppos <- grep(pattern=paste0("^Abundances Count [A-z0-9,. -]+",refchannel,"[A-z0-9,. -]+$"), names(data), value=FALSE)
  if (length(tmppos)) {
    collength <- length(names(data))
    if (PDversion>=21 & abdread) {
      data <- data[ ,c(1:(2*nread+5),tmppos)]
    } else {
      data <- data[ ,c(1:(nread+5),tmppos)]
    }
  } else {
    message("There is no protein abundance count information in the original data.")
    collength <- length(names(data))
    if (PDversion>=21 & abdread) {
      data <- data[ ,c(1:(2*nread+5))]
    } else {
      data <- data[ ,c(1:(nread+5))]
    }
    message("Set a default number for the protein abundance count to 2.")
    data$countNum <- 2
  }

  # Rename the channel to correct treatment
  for (i in 4:(nread+3)) {
    names(data)[i] <- treatment[i-3]
  }

  if (PDversion>=21 & abdread) {
    for (i in (nread+4):(2*nread+3)) {
      names(data)[i] <- paste0("Abundance_",treatment[i-(nread+3)])
    }
  }

  # Sum unique peptides, psms & Accession name correction
  names(data) <- gsub(pattern="^# Unique Peptides$","sumUniPeps", names(data))
  names(data) <- gsub(pattern="^Description$","description", names(data))
  names(data) <- gsub(pattern="^# PSMs$","sumPSMs", names(data))
  names(data) <- gsub(pattern=paste0("^Abundances Count [A-z0-9,. -]+",refchannel,"[A-z0-9,. -]+$"),"countNum", names(data))
  names(data) <- gsub(pattern="^Accession$","id", names(data))
  return(data)
}
