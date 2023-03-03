#' getName
#'
#' Functions to extract Gene/Protein/Both names from the Uniprot description
#'
#' @param x Uniprot description input
#' @param pfdatabase whether the input is from P.falciparum database
#'
#' @keywords internal
#'
#' @return character
#' @examples \dontrun{
#' }


getGeneName <- function(x, pfdatabase=FALSE) {
  if (pfdatabase) {
    gene=trimws(gsub("gene=", "", strsplit(x, "\\|")[[1]][2]))
  } else {
    gene=strsplit(strsplit(x, "GN=")[[1]][2], " ")[[1]][1]
  }
  if (length(gene)==0) {return(" ")} else {return(gene)}
}

getProteinName <- function(x, pfdatabase=FALSE) {
  if (pfdatabase) {
    protein=trimws(gsub("gene_product=", "", strsplit(x, "\\|")[[1]][4]))
  } else {
    protein=strsplit(x, " OS=")[[1]][1]
  }
  if (length(protein)==0) {return(" ")} else {return(protein)}
}

combineProteinGeneName <- function(data, printBothName, printGeneName, pfdatabase=FALSE) {
  if (printBothName) {
    data <- data %>% dplyr::rowwise() %>%
      dplyr::mutate(description1 = getProteinName(description, pfdatabase)) %>%
      dplyr::mutate(description2 = getGeneName(description, pfdatabase)) %>%
      dplyr::mutate(id = paste(id, description1, description2, sep="\n"))
    data$description1<-NULL
    data$description2<-NULL
  } else if (printGeneName) {
    data <- data %>% dplyr::rowwise() %>%
      dplyr::mutate(description = getGeneName(description, pfdatabase)) %>%
      dplyr::mutate(id = paste(id, description, sep="\n"))
  } else {
    data <- data %>% dplyr::rowwise() %>%
      dplyr::mutate(description = getProteinName(description, pfdatabase)) %>%
      dplyr::mutate(id = paste(id, description, sep="\n"))
  }
  data$description<-NULL
  return(data)
}
