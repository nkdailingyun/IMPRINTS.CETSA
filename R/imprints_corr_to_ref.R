#' imprints_corr_to_ref
#'
#' Function to find out the proteins with most similar profile to the reference profile
#'
#' @param data dataset of average IMPRINTS profile, that is after calculating the relative protein
#' abundance differences and deriving the average, i.e., after imprints_average(), make sure
#' the columns with readings are named in the format like "(Set_)37C_Treatment"
#' @param set a single character to indicate the sample name
#' @param treatment a single character to indicate the sample name
#' @param reference a numeric vector with the profile readings
#' @param use_score a single character element that define the method score. Method available : 'euclidean' or 'pearson'
#' @param score_threshold a numeric value to indicate the correlation score threshold, default set to 0.9
#' @param include_neg whether to include the negatively correlated proteins, default set to FALSE
#' @param max_na an integer indicating the maximum number of missing value for one protein, default is 0
#'
#'
#' @importFrom tidyr expand gather separate unite
#' @importFrom dplyr filter group_by left_join mutate rowwise summarise top_n ungroup
#' @importFrom tibble rownames_to_column
#' @export
#' @return a dataframe with the correlation from the profile and the proteins which passed the threshold
#' @examples \dontrun{
#' MOLM_M13_T1 <- imprints_corr_to_ref(MOLM, set="M13", treatment="T1",
#'   reference=c(0.02792811,0.03133724,0.14457743,0.33304728,0.27218208,0.23847792))
#' }
#'
#'

imprints_corr_to_ref <- function(data=NULL, set=NULL, treatment=NULL, reference=NULL,
                                 use_score = c("euclidean", "pearson"),
                                 score_threshold=0.9, include_neg=FALSE, max_na=0) {

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

  if(!(use_score %in% c("euclidean", "pearson"))){
    stop("Please choose a valid score method. For now, only 'euclidean' and 'pearson' are available.")
  }
  if (length(grep("countNum", names(data1)))) {
    countinfo1 <- unique(data1[ ,stringr::str_which(names(data1), "^id$|^sumPSM|^countNum|^sumUniPeps")])
    data1 <- data1[ ,-stringr::str_which(names(data1), "^sumPSM|^countNum|^sumUniPeps")]
  }
  if (length(grep("description", names(data1)))) {
    proteininfo <- unique(data1[ ,c("id","description")])
    data1$description <- NULL
  }
  nb_na <- apply(data1, 1, function(x) sum(is.na(x)))
  if (max_na < 0) {
    stop("Please enter an integer greater than or equal to 0 for the maximum NA per row")
  }
  na_thresh <- which(nb_na <= max_na)
  data1 <- data1[na_thresh, ]

  datal <- tidyr::gather(data1, treatment, reading, -id)
  if (length(unlist(strsplit(datal$treatment[1], "_")))==3) {
    datal <- tidyr::separate(datal, treatment, into=c("set","temperature","condition"))
    datal <- tidyr::unite(datal, "condition", c("set","condition"))
  }
  else if (length(unlist(strsplit(datal$treatment[1], "_")))==2) {
    datal <- tidyr::separate(datal, treatment, into=c("temperature","condition"))
  }

  treatment = unique(datal$condition)
  datal$condition <- NULL
  dataw <- tidyr::spread(datal, id, reading)
  dataw <- t(dataw)
  colnames(dataw) <- dataw[1, ]
  dataw <- dataw[-1, ]
  rname <- rownames(dataw)
  dataw <- apply(dataw, 2, as.numeric)
  rownames(dataw) <- rname
  dataw <- as.data.frame(dataw)

  if (use_score == "euclidean") { # introduced by Marc
    dataw$distance <- apply(dataw, 1, function(x) dist(rbind(x, reference))[1])
    dataw$score <- 1/(1 + dataw$distance)
    dataw$id <- rownames(dataw)
    rownames(dataw) <- 1:nrow(dataw)
    cortable <- dataw[ ,c("id", "score")]
    hist(cortable$score, breaks=100, xlab="", main="Distribution of Euclidean distance score")
  }
  else if (use_score == "pearson") {
    cortable <- cor(t(dataw[ ,1:ncol(dataw)]), reference, use = "pairwise.complete.obs")
    cortable <- as.data.frame(cortable)
    colnames(cortable) <- "score"
    cortable <- tibble::rownames_to_column(cortable, var = "id")
    hist(cortable$score, breaks=100, xlab="", main="Distribution of Pearson correlations")
  }

  if (score_threshold>0) {
    if (include_neg) {
      cortable <- subset(cortable, abs(correlation)>=score_threshold)
    } else {
      cortable <- subset(cortable, correlation>=score_threshold)
    }
  }
  if (nrow(cortable)==0) {
    # stop(paste0("No proteins pass the correlation threshold of ", score_threshold, "\n try to lower the threshold..."))
    g <- ggplot(data.frame(x = c(0,1), y = c(0,1)), aes(x,y, label = "s")) +
      geom_text(x=0.5, y=0.5, label = paste0("No proteins pass the correlation score threshold of ",
                                             score_threshold, "\n try to lower the threshold",
                                             "\nor change the maximum number",
                                             "\nof missing values"), size = 6) +
      cowplot::theme_cowplot() +
      theme(axis.text.x = element_blank(),
            axis.title.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.text.y = element_blank(),
            axis.title.y = element_blank(),
            axis.ticks.y = element_blank())

    return(g)
  } else {
    message(paste0(nrow(cortable), " proteins pass the correlation score threshold of ", score_threshold))
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
  return(cortable[order(cortable$score,decreasing=TRUE), ])
}
