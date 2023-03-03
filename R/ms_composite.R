#' ms_composite_ID_Gene
#'
#' Functions to make a composite Uniprot ID/Gene name from the Uniprot description
#'
#' @param data a dataframe with id and description column
#' @param pfdatabase whether the input is from P.falciparum database
#'
#' @export
#'
#' @return a dataframe with a compositie id column
#' @examples \dontrun{
#' }

ms_composite_ID_Gene <- function(data, pfdatabase=FALSE) {
  data <- data %>% dplyr::rowwise() %>%
    dplyr::mutate(description = getGeneName(description, pfdatabase)) %>%
    dplyr::mutate(id = paste(id, description, sep="\n"))
  data$description<-NULL
  return(data)
}

#' ms_composite_ID_Protein
#'
#' Functions to make a composite Uniprot ID/Protein name from the Uniprot description
#'
#' @param data a dataframe with id and description column
#' @param pfdatabase whether the input is from P.falciparum database
#'
#' @export
#'
#' @return a dataframe with a compositie id column
#' @examples \dontrun{
#' }

ms_composite_ID_Protein <- function(data, pfdatabase=FALSE) {
  data <- data %>% dplyr::rowwise() %>%
    dplyr::mutate(description = getProteinName(description, pfdatabase)) %>%
    dplyr::mutate(id = paste(id, description, sep="\n"))
  data$description<-NULL
  return(data)
}

#' ms_composite_ID_Gene_Protein
#'
#' Functions to make a composite Uniprot ID/Gene name/Protein name from the Uniprot description
#'
#' @param data a dataframe with id and description column
#' @param pfdatabase whether the input is from P.falciparum database
#'
#' @export
#'
#' @return a dataframe with a compositie id column
#' @examples \dontrun{
#' }
#'

ms_composite_ID_Gene_Protein <- function(data, pfdatabase=FALSE) {
  data <- data %>% dplyr::rowwise() %>%
    dplyr::mutate(description1 = getProteinName(description, pfdatabase)) %>%
    dplyr::mutate(description2 = getGeneName(description, pfdatabase)) %>%
    dplyr::mutate(id = paste(id, description1, description2, sep="\n"))
  data$description1<-NULL
  data$description2<-NULL
  data$description<-NULL
  return(data)
}

#' ms_getGene
#'
#' Functions to add one extra column with gene information from the Uniprot description
#'
#' @param data a dataframe with id and description column
#' @param pfdatabase whether the input is from P.falciparum database
#'
#' @export
#'
#' @return a dataframe with an extra gene column
#' @examples \dontrun{
#' }
#'
ms_getGene <- function(data, pfdatabase=FALSE) {
  data <- data %>% dplyr::rowwise() %>%
    dplyr::mutate(gene = getGeneName(description, pfdatabase))
  return(data)
}
