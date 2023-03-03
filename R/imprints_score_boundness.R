#' imprints_score_boundness
#'
#' Function to categorize protein according to their Expression/Stability change
#'
#' @param data dataset after calculating the relative protein abundance differences
#' i.e., imprints_caldiff()
#' @param format choose between 9 and 4, indicating how many categories to segregate,
#' default value is 9, in which proteins are divided into 9 categories: NN, CN+, CN-,
#' NC+, NC-, CC++, CC+-, CC-+, CC--. When switch to 4, proteins are instead divided into
#' 4 categories: NN, CN, NC, CC. The sign of + or - after N or C is determined by the
#' sign of abundance.score and stability.score.mean, respectively
#' @param basetemp character indicating the baseline temperature for expression levels, default value is 37C
#' @param meancutoff The threshold value for the mean of triplicate measurement (First filter of hits)
#' @param boundedness Number of SEM, used for setting a bound that the mean needs to be within
#' @param qualitycutoff The threshold value for the SEM to separate well-measured and ill-measured measurements
#' @param meandev The deviation of individual mean from the group mean in CN categorization
#'
#' @import stringr
#' @export
#' @return a list
#' @examples \dontrun{
#'   MOLM <- imprints_score_boundness(MOLM)
#' }
#'
#'


imprints_score_boundness <- function(data, format="9", basetemp="37C", meancutoff=0.25,
                                     boundedness=4, qualitycutoff=0.15, meandev=0.15) {

  dataname <- deparse(substitute(data))
  outdir <- ms_directory(data, dataname)$outdir
  data <- ms_directory(data, dataname)$data

  if (!as.character(format) %in% c("4","9")) {
    stop("Please specify format=9 or format=4 !")
  }

  if (length(grep("description", names(data)))) {
    proteininfo <- unique(data[ ,c("id","description")])
    data$description <- NULL
  }
  if (length(grep("countNum", names(data)))) {
    countinfo <- unique(data[ ,stringr::str_which(names(data), "^id$|^sumPSM|^countNum|^sumUniPeps")])
    data <- data[ ,-stringr::str_which(names(data), "^sumPSM|^countNum|^sumUniPeps")]
  }

  refcol <- which(apply(data[,-1],2,sum,na.rm=T)==0)
  data <- tidyr::gather(data[ ,-(refcol+1)], condition, reading, -id)
  # data <- gather(data, condition, reading, -id)
  if (length(unlist(strsplit(data$condition[1], "_")))==4) {
    data <- tidyr::separate(data, condition, into=c("set","temperature","replicate","treatment"), sep="_")
    data_clean <- data %>%
      dplyr::group_by(id,set,temperature,treatment) %>%
      dplyr::summarise(Mean = mean(reading,na.rm=T),
                      SD = sd(reading,na.rm=T),
                      SEM = SD/sqrt(length(na.omit(reading))))

    data_abd <- data %>% dplyr::group_by(id, set, treatment) %>%
      dplyr::filter(temperature==basetemp) %>%
      dplyr::summarise(abundance.score=mean(reading,na.rm=T))
    data_thermal <- data %>% dplyr::left_join(data_abd) %>% dplyr::rowwise() %>%
      dplyr::mutate(stability=reading-abundance.score) %>% dplyr::filter(temperature!=basetemp)
    data_thermal1 <- data_thermal %>% dplyr::group_by(id,set,treatment,temperature) %>%
      dplyr::summarise(stability.score=mean(stability,na.rm=T))
    data_thermal2 <- data_thermal1 %>% dplyr::group_by(id,set,treatment) %>%
      dplyr::summarise(stability.score.mean=mean(stability.score,na.rm=T))
  } else if (length(unlist(strsplit(data$condition[1], "_")))==3) {
    data <- tidyr::separate(data, condition, into=c("temperature","replicate","treatment"), sep="_")
    data_clean <- data %>%
      dplyr::group_by(id,temperature,treatment) %>%
      dplyr::summarise(Mean = mean(reading,na.rm=T),
                      SD = sd(reading,na.rm=T),
                      SEM = SD/sqrt(length(na.omit(reading))))

    data_abd <- data %>% dplyr::group_by(id, treatment) %>%
      dplyr::filter(temperature==basetemp) %>%
      dplyr::summarise(abundance.score=mean(reading,na.rm=T))
    data_thermal <- data %>% dplyr::left_join(data_abd) %>% dplyr::rowwise() %>%
      dplyr::mutate(stability=reading-abundance.score) %>% dplyr::filter(temperature!=basetemp)
    data_thermal1 <- data_thermal %>% dplyr::group_by(id,treatment,temperature) %>%
      dplyr::summarise(stability.score=mean(stability,na.rm=T))
    data_thermal2 <- data_thermal1 %>% dplyr::group_by(id,treatment) %>%
      dplyr::summarise(stability.score.mean=mean(stability.score,na.rm=T))
  } else {
    stop("make sure the namings of the columns of the dasaset are correct.")
  }

  data_thermal3 <- merge(data_abd, data_thermal2)
  data_thermal3 <- data_thermal3 %>% dplyr::rowwise() %>%
    dplyr::mutate(abundance.sign=ifelse(abundance.score>0,T,F))
  data_thermal3 <- data_thermal3 %>% dplyr::rowwise() %>%
    dplyr::mutate(stability.sign=ifelse(stability.score.mean>0,T,F))
  data_thermal3 <- data_thermal3 %>%
    dplyr::mutate(category=paste0(as.character(abundance.sign),as.character(stability.sign)))
  data_thermal3$category <- gsub("FALSE","N",data_thermal3$category)
  data_thermal3$category <- gsub("TRUE","C",data_thermal3$category)
  data_thermal3 <- data_thermal3[ ,c("id","treatment","category")] %>% ungroup()

  selection_metrics <- data_clean %>%
    dplyr::group_by(mean_threshold = (abs(Mean) > meancutoff),
                    bounded = (abs(Mean)- boundedness*SEM>0),
                    wellmeasured = (SEM < qualitycutoff))

  ## Initial hitlist (Reference)
  # Definition for hit
  hits_definition <- selection_metrics %>%
    dplyr::filter(mean_threshold == T, bounded == T)

  #Unique id, treatment pair for hits
  if (length(grep("set",names(data_thermal2)))==0) {
    ##### need to implement the one with set possibility in future ######
    hits_keys <- hits_definition %>% dplyr::ungroup() %>%
      dplyr::select(id,treatment) %>% dplyr::distinct()

    # Reference hitlist
    hitlist <- selection_metrics %>% dplyr::right_join(hits_keys, by=names(hits_keys))
    referencelist <- hitlist

    # Separate hits and NN from data
    NN <- selection_metrics %>% dplyr::anti_join(hitlist) %>% dplyr::group_by(category = 'NN')

    ## Extract Not determinable -  ND
    # Proteins with noisy at baseline temperature
    ND_keys <- hitlist %>% dplyr::filter(temperature == basetemp,
                                         bounded == FALSE,
                                         wellmeasured == FALSE) %>%
               dplyr::ungroup() %>% dplyr::select(id,treatment) %>% dplyr::distinct()
    ND <- hitlist %>% dplyr::right_join(ND_keys, by=names(ND_keys))

    # Proteins without measurement at baseline temperature, combined into ND
    NDwo37 <- hitlist %>% dplyr::filter(temperature == basetemp,
                                        is.na(Mean) == TRUE)
    ND <- ND %>% dplyr::full_join(NDwo37, by=names(ND)) %>% dplyr::group_by(category = 'ND')
    ND <- dplyr::inner_join(ND, data_thermal3, by=c("id","treatment")) %>%
          dplyr::mutate(category=paste0(category.x,category.y)) %>%
          dplyr::select(-c(category.x,category.y))

    # Remove ND from hitlist
    hitlist <- hitlist %>% dplyr::anti_join(ND, by=c('id','treatment'))

    ## Stability change w/out expression change  - NC
    # High Temp. hits
    NC_highT <- hitlist %>% dplyr::filter(temperature != basetemp,
                                          mean_threshold == TRUE,
                                          bounded == TRUE) %>%
      dplyr::ungroup() %>% dplyr::select(id, treatment) %>% dplyr::distinct()

    # Low Temp., small mean and well measured
    NC_lowT <- hitlist %>% dplyr::filter(temperature == basetemp,
                                         mean_threshold == FALSE,
                                         wellmeasured == TRUE) %>%
      dplyr::ungroup() %>% dplyr::select(id,treatment) %>% distinct()

    # Match id,condition pair fulfilling NC condition
    NC_keys <- dplyr::inner_join(NC_highT, NC_lowT)
    NC <- hitlist %>% dplyr::right_join(NC_keys, names(NC_keys)) %>% dplyr::group_by(category='NC')
    NC <- NC %>% dplyr::inner_join(data_thermal3, by=c("id","treatment")) %>%
      dplyr::mutate(category=paste0(category.x,category.y)) %>%
      dplyr::select(-c(category.x,category.y))

    # Remove NC from hitlist
    hitlist <- hitlist %>% dplyr::anti_join(NC, by=c('id','treatment'))
  }

  ## Expression change w/out stability change - CN

  ################################## INTERVAL CONDITION FOR CN ##################################
  CN_intervals <- hitlist %>% dplyr::group_by(upperbound=Mean+SEM, lowerbound=Mean-SEM)
  CN_keys <- CN_intervals %>% dplyr::ungroup() %>% dplyr::select(id, treatment) %>% dplyr::distinct()

  CN_keep_interval <- c()
  CN_discard_interval <- c() # Will be part of CC

  for (idx in 1:nrow(CN_keys)) {
    iterTibble <- CN_intervals %>%
      dplyr::filter(id == CN_keys$id[idx], treatment == CN_keys$treatment[idx]) %>%
      dplyr::select(temperature,upperbound,lowerbound)
    overlap <- c()
    # print(iterTibble)
    for (i in 1:nrow(iterTibble)) {
      refinterval <- iterTibble[i, ]
      testinterval <- iterTibble[-i, ]
      overlap <- c(overlap,sum(refinterval$lowerbound <= testinterval$upperbound & refinterval$upperbound >= testinterval$lowerbound)>0)
    }
    if (all(overlap, na.rm=T)) {
      CN_keep_interval <- rbind(CN_keep_interval,CN_keys[idx,])
    } else {
      CN_discard_interval <- rbind(CN_discard_interval,CN_keys[idx,])
    }
  }
  ####################################################################################################

  ############################## DEVIATION FROM REFERENCE MEAN CONDITION #############################
  # baseline temperature measurement with high mean
  CN_referencekeys <- hitlist %>%
    dplyr::filter(temperature == basetemp, mean_threshold == TRUE) %>%
    dplyr::ungroup () %>% dplyr::select(id,treatment) %>% dplyr::distinct()

  CN_keep_meandev <- c()
  CN_discard_meandev <- c() # Will be part of CC

  for (idx in 1:nrow(CN_referencekeys)) {
    iterTibble <- hitlist %>% dplyr::filter(id == CN_referencekeys$id[idx],
                                     treatment == CN_referencekeys$treatment[idx])
    if (nrow(iterTibble %>% dplyr::filter(temperature != basetemp, bounded == TRUE)) == 0) {
      next # Not bounded high temps are disregarded
    } else {
      refmean <- (iterTibble %>% dplyr::filter(temperature == basetemp))$Mean
      boundediterTibble <- iterTibble %>% dplyr::filter(temperature != basetemp, bounded == TRUE)
      deviateTibble <- iterTibble %>% dplyr::filter(temperature != basetemp, wellmeasured == TRUE)
      # Bounded high temp. deviate at most meandev AND
      # wellmeasured high temp don't deviate more than meandev
      if ( (sum((abs(boundediterTibble$Mean - refmean) < meandev*abs(refmean))) == nrow(boundediterTibble)) &&
           (!any(abs(deviateTibble$Mean - refmean) > meandev*abs(refmean))) ) {
        CN_keep_meandev <- rbind(CN_keep_meandev,CN_referencekeys[idx,])
      } else {
        CN_discard_meandev <- rbind(CN_discard_meandev,CN_referencekeys[idx,])
      }
    }
  }
  ######################################################################################################

  # Proteins that do not fulfill either of CN conditions are CC - Expression and stability change
  CC_keys <- dplyr::inner_join(CN_discard_interval, CN_discard_meandev)
  CC <- hitlist %>% dplyr::right_join(CC_keys,by=names(CC_keys)) %>% dplyr::group_by(category = 'CC')
  CC <- CC %>% dplyr::inner_join(data_thermal3, by=c("id","treatment")) %>%
    dplyr::mutate(category=paste0(category.x,category.y)) %>%
    dplyr::select(-c(category.x,category.y))

  # Remove all CC and reveal CNs
  CN <- hitlist %>% dplyr::anti_join(CC, by=c('id','treatment')) %>% dplyr::group_by(category = 'CN')
  CN <- CN %>% dplyr::inner_join(data_thermal3, by=c("id","treatment")) %>%
    dplyr::mutate(category=paste0(category.x,category.y)) %>%
    dplyr::select(-c(category.x,category.y))

  ##################################################
  #                                                #
  #    Data Wrangling: Compiling and exporting     #
  #                                                #
  ##################################################

  # Complete set of all categorized proteins
  categorized_full <- CC %>%
    dplyr::full_join(CN, by = colnames(CN)) %>%
    dplyr::full_join(NC, by = colnames(NC)) %>%
    dplyr::full_join(ND, by = colnames(ND)) %>% distinct()

  if (as.character(format)=="9") {
    categorized_full$category <- gsub("ND[CN]{2}","ND",categorized_full$category)
    categorized_full$category <- gsub("NN[CN]{2}","NN",categorized_full$category)
    categorized_full$category <- gsub("CNC[CN]{1}","CN+",categorized_full$category)
    categorized_full$category <- gsub("CNN[CN]{1}","CN-",categorized_full$category)
    categorized_full$category <- gsub("NC[CN]{1}N","NC-",categorized_full$category)
    categorized_full$category <- gsub("NC[CN]{1}C","NC+",categorized_full$category)
    categorized_full$category <- gsub("CCCC","CC++",categorized_full$category)
    categorized_full$category <- gsub("CCCN","CC+-",categorized_full$category)
    categorized_full$category <- gsub("CCNC","CC-+",categorized_full$category)
    categorized_full$category <- gsub("CCNN","CC--",categorized_full$category)
  }
  else if (as.character(format)=="4") {
    categorized_full$category <- gsub("ND[CN]{2}","ND",categorized_full$category)
    categorized_full$category <- gsub("NN[CN]{2}","NN",categorized_full$category)
    categorized_full$category <- gsub("CN[CN]{2}","CN",categorized_full$category)
    categorized_full$category <- gsub("NC[CN]{2}","NC",categorized_full$category)
    categorized_full$category <- gsub("CC[CN]{2}","CC",categorized_full$category)
  }

  out_reference <- dplyr::inner_join(proteininfo, referencelist, by = 'id')
  out_categorized <- dplyr::inner_join(proteininfo, categorized_full, by = 'id')
  out_summary <- categorized_full %>% dplyr::group_by(id, treatment, category) %>% dplyr::summarise()
  out_summary <- dplyr::inner_join(proteininfo, out_summary, by = 'id')
  out_NN <- dplyr::inner_join(proteininfo, NN, by = 'id')

  ms_filewrite(out_reference, paste0(dataname, "_boundness_Hitlist.txt"), outdir=outdir)
  ms_filewrite(out_categorized, paste0(dataname, "_boundness_Categorized.txt"), outdir=outdir)
  ms_filewrite(out_summary, paste0(dataname, "_boundness_Summary.txt"), outdir=outdir)
  ms_filewrite(out_NN, paste0(dataname, "_boundness_NN.txt"), outdir=outdir)

  ##################################################
  #                                                #
  #          Summary: Return to user               #
  #                                                #
  ##################################################

  ### Return : Summary of run and data frame with categories for further handling. ###
  summary_NN <- NN %>% dplyr::group_by(id,treatment,category) %>% dplyr::summarise()
  tablecate <- table(out_summary[ ,c(3:4)])
  tableNN <- table(summary_NN[ ,c(2:3)])
  message('This dataset has the following count for the categories:')
  message(cbind(tablecate,tableNN))

  message('The categories were obtained with the following values for the parameters:')
  message(sprintf(fmt = 'meancutoff: %s', meancutoff))
  message(sprintf(fmt = 'boundedness: %s', boundedness))
  message(sprintf(fmt = 'qualitycutoff: %s', qualitycutoff))
  message(sprintf(fmt = 'meandev: %s', meandev))

  results <- rbind(out_categorized, out_NN)
  results.list <- split(results, f = results$category)
  return(results.list)
}
