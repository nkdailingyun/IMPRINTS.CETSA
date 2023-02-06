#' imprints_corr_to_ref
#'
#' Function to find out the proteins with most similar profile to the reference profile,
#' using the Pearson correlation metric
#'
#' @param data dataset of average IMPRINTS profile, that is after calculating the relative protein
#' abundance differences and deriving the average, i.e., after imprints_average(), make sure
#' the columns with readings are named in the format like "(Set_)37C_Treatment"
#' @param set a single character to indicate the sample name
#' @param treatment a single character to indicate the sample name
#' @param reference a numeric vector with the profile readings
#' @param corr_threshold a numeric value to indicate the threshold, default set to 0.8
#' @param include_neg whether to include the negatively correlated proteins,
#' default set to FALSE
#'
#'
#' @import dplyr
#' @export
#' @return a dataframe
#' @examples \dontrun{
#' MOLM_M13_T1 <- imprints_corr_to_ref(MOLM, set="M13", treatment="T1",
#'   reference=c(0.02792811,0.03133724,0.14457743,0.33304728,0.27218208,0.23847792))
#' }
#'
#'

imprints_corr_to_ref <- function(data=NULL, set=NULL, treatment=NULL, reference=NULL,
                              corr_threshold=0.8, include_neg=FALSE) {

  if (!length(set) %in% c(0,1)) {stop("Pls provide one set name if any")}
  if (length(treatment) != 1) {stop("Pls provide one treatment name")}
  if (length(reference) == 0) {stop("Pls provide a numeric reading vector as the reference profile")}

  if (inherits(data,"data.frame") | class(data)=="data.frame") {
    dataname <- deparse(substitute(data))
    outdir <- ms_directory(data, dataname)$outdir
    data <- ms_directory(data, dataname)$data

    subset <- grep(treatment,names(data))
    if (length(set)) {
      subset_set <- grep(set,names(data))
      subset <- intersect(subset,subset_set)
    }
    if (length(subset)>0 & length(subset)==length(reference)) {
      data1 = data[ ,c(1,2,subset)]
    } else {
      stop("Pls provide the right set/treatment keyword character")
    }
  }

  if (length(grep("countNum", names(data1)))) {
    countinfo1 <- unique(data1[ ,c("id","sumUniPeps","sumPSMs","countNum")])
    data1 <- data1[ ,!(names(data1) %in% c("sumUniPeps","sumPSMs","countNum"))]
  }
  if (length(grep("description", names(data1)))) {
    proteininfo <- unique(data1[ ,c("id","description")])
    data1$description <- NULL
  }

  datal <- gather(data1, treatment, reading, -id)
  if (length(unlist(strsplit(datal$treatment[1], "_")))==3) {
    datal <- separate(datal, treatment, into=c("set","temperature","condition"))
    datal <- unite(datal, "condition", c("set","condition"))
  } else if (length(unlist(strsplit(datal$treatment[1], "_")))==2) {
    datal <- separate(datal, treatment, into=c("temperature","condition"))
  }
  treatment = unique(datal$condition)
  datal$condition <- NULL
  dataw <- spread(datal, id, reading)
  # return(datal)
  # return(dataw)
  #print(as.character(eval(parse(text = treatment[1]))))
  cortable <- cor(as.matrix(dataw[ ,-1]), reference, use="pairwise.complete.obs")
  cortable <- as.data.frame(cortable)
  colnames(cortable) <- "correlation"
  cortable <- tibble::rownames_to_column(cortable, var = "id")
  # return(cortable)
  hist(cortable$correlation, breaks=100, xlab="", main="Distribution of Pearson correlations")
  if (corr_threshold>0) {
    if (include_neg) {
      cortable <- subset(cortable, abs(correlation)>=corr_threshold)
    } else {
      cortable <- subset(cortable, correlation>=corr_threshold)
    }
  }
  if (nrow(cortable)==0) {
    stop(paste0("No proteins pass the correlation threshold of ", corr_threshold, "\n try to lower the threshold..."))
  } else {
    print(paste0(nrow(cortable), " proteins pass the correlation threshold of ", corr_threshold))
  }
  # cortable <- group_by(dataw,gene) %>%
  #   summarize(correlation=cor(eval(parse(text=treatment[1])),eval(parse(text=treatment[2])), use="complete.obs"))
  # dataw <- spread(datal, temperature, reading)
  # dism <- plyr::ddply(dataw, "gene", function(data) {
  #   data<-data[order(data$condition), ]
  #   dm<-as.matrix(dist(data[ ,-c(1:2)]))[1,2]
  # })
  # names(dism)[2] <- "distance"
  cortable <- merge(proteininfo,cortable)
  # cortable <- merge(cortable,dism)
  return(cortable[order(cortable$correlation,decreasing=T), ])
}
