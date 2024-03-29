#' imprints_rawread
#'
#' Function to parse and read in IMPRINTS-CETSA data from tab delimited files exported from Proteome Discoverer
#'
#' @param filevector a file name or a vector of filenames to import
#' @param fchoose whether to choose file interactively, default set to FALSE
#' @param treatment a vector of treatment names applied to CETSA samples, in the same order as channels,
#' in a typical IMPRINTS/2D working scheme, the treatment name should contain both replicate and treatment
#' information, preferably concatenated with an underline
#' @param nread number of reading channels, should match the number of channels used, default value 10
#' @param fdrcontrol whether to check the protein FDR confidence level, default set to FALSE
#' @param refchannel names of reference channel used in Proteome Discoverer search, default value 126
#' @param channels names of the read-in channels, default value NULL, it would automatically
#' match the provided channel number when it is 10, 11, 16 or 18
#' @param software software used to process the IMPRINTS-CETSA raw data and obtain ProteinGroup files.
#' Either Proteome Discoverer (PD) or MaxQuant is currently supported. Default is PD.
#'
#' @importFrom dplyr filter group_by left_join mutate rowwise summarise top_n ungroup
#' @importFrom magrittr %>%
#' @export
#' @return a dataframe
#' @examples \dontrun{
#'  MOLM13 <- imprints_rawread(c("M13_37C_frac_Proteins.txt","M13_52C_frac_Proteins.txt","M13_58C_frac_Proteins.txt"), nread=10,
#'  treatment=("B1_DMSO","B1_TNFa","B1_AT26533","B2_DMSO","B2_TNFa","B2_AT26533","B3_DMSO","B3_TNFa","B3_AT26533","Mix"))
#' }
#'
#'
#'
imprints_rawread <- function(filevector, fchoose=FALSE, treatment=NULL, nread=10,
                             fdrcontrol=FALSE, refchannel="126", channels=NULL,
                             software = c("PD", "MaxQuant", "pFind")) {

  software <- match.arg(software)
  if (nread==10 & length(channels)==0) {
    channels=c("126","127N","127C","128N","128C","129N","129C","130N","130C","131")
  } else if (nread==11 & length(channels)==0) {
    channels=c("126","127N","127C","128N","128C","129N","129C","130N","130C","131N","131C")
  } else if (nread==16 & length(channels)==0) {
    channels=c("126","127N","127C","128N","128C","129N","129C","130N","130C","131N","131C",
               "132N","132C","133N","133C","134N")
  } else if (nread==18 & length(channels)==0) {
    channels=c("126","127N","127C","128N","128C","129N","129C","130N","130C","131N","131C",
               "132N","132C","133N","133C","134N","134C","135N")
  } else if (nread!=length(channels) | nread!=length(treatment)) {
    stop("Please provide a vector of used TMT channels")
  }

  if (length(treatment)==0) {
    stop("Need to specify the treatment conditions in the same order as the TMT channel arrangement")
  }
  if (length(treatment)!=nread | length(treatment)!=length(channels)) {
    stop("The number of elements in treatment condition vector should be equal to the number of read-in channels")
  }
  if (length(filevector)==0) {
    stop("Need to specify the input data file names")
  }
  flength <- length(filevector)

  if (flength < 2) {
    dirname <- deparse(substitute(filevector))
    dirname_l <- unlist(strsplit(dirname, split="/"))
    dirname <- dirname_l[length(dirname_l)]
    data <- ms_innerread(filevector, fchoose, treatment, nread,
                         fdrcontrol, refchannel, channels, software)
    data <- ms_dircreate(dirname, data)
    outdir <- attr(data,"outdir")
    if (length(attr(data,"outdir"))==0 & length(outdir)>0) {
      attr(data,"outdir") <- outdir
    }
    message("The data composition under each experimental condition (read in) is:")
    print(table(data$condition))
    return(data)
  }
  else {
    filename <- filevector[1]
    dirname <- deparse(substitute(filename))
    dirname_l <- unlist(strsplit(dirname, split="/"))
    dirname <- dirname_l[length(dirname_l)]
    indata <- ms_innerread(filevector[1], fchoose, treatment, nread,
                           fdrcontrol, refchannel, channels, software)
    indata <- dplyr::mutate(indata, condition = paste0(condition,".1"))
    outdata <- indata
    for (i in 2:flength) {
      indata <- ms_innerread(filevector[i], fchoose, treatment, nread,
                             fdrcontrol, refchannel, channels, software)
      indata <- dplyr::mutate(indata, condition = paste0(condition, ".", i))
      outdata <- rbind(x=outdata, y=indata)
    }
    outdata <- ms_dircreate(paste0("merged_",dirname), outdata)
    outdir <- attr(outdata,"outdir")
    if (length(attr(outdata,"outdir"))==0 & length(outdir)>0) {
      attr(outdata,"outdir") <- outdir
    }

    message("The data composition under each experimental condition (read in) is:")
    print(table(outdata$condition))
    rownames(outdata) <- NULL

    return(outdata)
  }
}
