#' imprints_normalization
#'
#' Function to perform a systematic normalization on the whole dataset
#' First IMPRINTS-CETSA data is transformed into an eSet format, then a VSN-based model normalization is applied,
#' and finally the function ‘removeBatchEffect’ from limma R package is called to remove potential batch effect.
#'
#' @param datafile a directory pathway pointing to a file with pre-normalization
#' data and experimental scheme information, under the following keywords:
#' set, replicate, temperature, treatment
#' @param isfile whether the input datafile is a file, which should contain the
#' annotation rows on the top of each quantitative column of the saved dataset
#' with keyword information of “temperature”, “treatment”, “replicate”, or
#' “set” in the case of complex experimental layout, default set to FALSE
#' @param todf whether to return the normalized data in a dataframe format,
#' default set to TRUE
#' @param pertemp whether do the normalization on each heating temperature separately,
#' default set to TRUE
#'
#' @import vsn
#' @importFrom limma contrasts.fit eBayes lmFit makeContrasts topTable
#' @importFrom Biobase exprs phenoData
#' @export
#' @return a dataframe
#' @examples \dontrun{
#'  MOLM <- imprints_normalization(MOLM)
#' }
#'
#'

imprints_normalization <- function(datafile, isfile=FALSE, todf=TRUE, pertemp=TRUE) {

  if (isfile) {
    outdir <- stringr:str_split(datafile,"/")[[1]][2]
    rows <- read.table(file=datafile, header=FALSE, nrows=10, fill=T)[ ,1]
    pos <- grep("^id$", rows)
    data <- readr::read_tsv(file=datafile, skip=(pos-1))
    pdata <- NULL
    pdata <- readr::read.table(file=datafile, header=FALSE, nrows=(pos-1), stringsAsFactors=FALSE, row.names=NULL)
    if (nrow(pdata)==0) {
      stop("Remember to add the necessary information of experimental scheme in the data file!")
    }
    nread <- ncol(pdata)-1
    pdata1 <- data.frame(t(pdata[ ,-1]), row.names=names(data[ ,c(3:(nread+2))]))
    names(pdata1) <- pdata[ ,1]
    #return(list(data=data, metadata=metadata1))
    data_eset <- ms_to_eSet(data=data, nread=nread, pdata=pdata1)
  } else {
    data <- datafile
    data <- data[order(data$id), ] #important to first sort the IDs
    dataname <- deparse(substitute(datafile))
    outdir <- ms_directory(datafile, dataname)$outdir
    data <- ms_directory(datafile, dataname)$data
    cname <- setdiff(names(datafile), c("id","description","sumUniPeps","sumPSMs","countNum"))
    if (length(unlist(strsplit(cname[1], "_")))==3) {
      temperature <- unlist(lapply(strsplit(cname, "_"),`[`,1))
      replicate <- unlist(lapply(strsplit(cname, "_"),`[`,2))
      treatment <- unlist(lapply(strsplit(cname, "_"),`[`,3))
      pdata1 <- data.frame(temperature=temperature, replicate=replicate, treatment=treatment)
      row.names(pdata1) <- cname
      nread <- nrow(pdata1)
      data_eset <- ms_to_eSet(data=data, nread=nread, pdata=pdata1)
    } else if (length(unlist(strsplit(cname[1], "_")))==4) {
      set <- unlist(lapply(strsplit(cname, "_"),`[`,1))
      temperature <- unlist(lapply(strsplit(cname, "_"),`[`,2))
      replicate <- unlist(lapply(strsplit(cname, "_"),`[`,3))
      treatment <- unlist(lapply(strsplit(cname, "_"),`[`,4))
      pdata1 <- data.frame(set=set, temperature=temperature, replicate=replicate, treatment=treatment)
      row.names(pdata1) <- cname
      nread <- nrow(pdata1)
      data_eset <- ms_to_eSet(data=data, nread=nread, pdata=pdata1)
    } else {
      stop("make sure the namings of the columns of the dasaset are correct.")
    }
  }

  data_eset_mock <- ms_to_eSet(data=na.omit(data), nread=nread, pdata=pdata1)

  try(arrayQualityMetrics::arrayQualityMetrics(expressionset=data_eset_mock,
                                               outdir=paste0(dataname,"_pre_norm_QC_report"),
                                               force=TRUE, do.logtransform=TRUE,
                                               intgroup=Biobase::phenoData(data_eset_mock)@varMetadata$labelDescription,
                                               reporttitle=paste0(dataname, " Pre Normalization QC report")))

  # to perform vsn normalization, whether per temperature can be controlled
  if (pertemp) {
    data1 <- Biobase::exprs(data_eset)
    if (sum(grepl("set",names(pdata1)))) {
      uniset <- unique(as.character(data_eset$set))
      unitemp <- unique(as.character(data_eset$temperature))
      for (se in uniset) {
        for (temp in unitemp) {
          pos <- intersect(grep(se, colnames(data1)),grep(temp, colnames(data1)))
          data_sub <- Biobase::exprs(data_eset[ ,data_eset$set==se & data_eset$temperature==temp])
          if (nrow(data_sub)<43) {
            data_sub <- vsn::justvsn(data_sub,
                                     minDataPointsPerStratum=nrow(data_sub)-1)
          } else {
            data_sub <- vsn::justvsn(data_sub)
          }
          data1[ ,pos] <- data_sub
        }
      }
      # print(head(data1))
    } else {
      unitemp <- unique(as.character(data_eset$temperature))
      for (temp in unitemp) {
        pos <- grep(temp, colnames(data1))
        data_sub <- Biobase::exprs(data_eset[ ,data_eset$temperature==temp])
        if (nrow(data_sub)<43) {
          data_sub <- vsn::justvsn(data_sub,
                                   minDataPointsPerStratum=nrow(data_sub)-1)
        } else {
          data_sub <- vsn::justvsn(data_sub)
        }
        data1[ ,pos] <- data_sub
      }
      # print(head(data1))
    }
    data_eset <- ms_to_eSet(data1, matrixonly=TRUE,
                            nread=ncol(data1), refeset=data_eset)
  } else { # do not do per temperature normalization
    if (sum(grepl("set",names(pdata1)))) {
      data1 <- Biobase::exprs(data_eset)
      uniset <- unique(as.character(data_eset$set))
      for (se in uniset) {
        pos <- grep(se, colnames(data1))
        data_sub <- Biobase::exprs(data_eset[ ,data_eset$set==se])
        if (nrow(data_sub)<43) {
          data_sub <- vsn::justvsn(data_sub,
                                   minDataPointsPerStratum=nrow(data_sub)-1)
        } else {
          data_sub <- vsn::justvsn(data_sub)
        }
        data1[ ,pos] <- data_sub
      }
      data_eset <- ms_to_eSet(data1, matrixonly=TRUE,
                              nread=ncol(data1), refeset=data_eset)
    } else {
      if (nrow(Biobase::exprs(data_eset))<43) {
        data_eset <- vsn::justvsn(data_eset,
                                  minDataPointsPerStratum=nrow(Biobase::exprs(data_eset))-1)
      } else {
        data_eset <- vsn::justvsn(data_eset)
      }
    }
  }
  data_eset1 <- limma::removeBatchEffect(Biobase::exprs(data_eset),
                                         batch=as.character(pData(data_eset)$replicate),
                                         design=model.matrix(~pData(data_eset)$treatment))
  # to perform quantile normalization
  # if (sum(grepl("set",names(pdata1)))) {
  #   for (se in uniset) {
  #     for (temp in unitemp) {
  #       pos <- intersect(grep(se, colnames(data_eset1)),grep(temp, colnames(data_eset1)))
  #       data_sub <- data_eset1[ ,pos]
  #       data_sub <- limma::normalizeQuantiles(data_sub)
  #       data_eset1[ ,pos] <- data_sub
  #     }
  #   }
  #   print(head(data_eset1))
  # } else {
  #   for (temp in unitemp) {
  #     pos <- grep(temp, colnames(data_eset1))
  #     data_sub <- data_eset1[ ,pos]
  #     data_sub <- limma::normalizeQuantiles(data_sub)
  #     data_eset1[ ,pos] <- data_sub
  #   }
  #   print(head(data_eset1))
  # }

  data_eset1 <- ms_to_eSet(data_eset1, matrixonly=TRUE,
                           nread=ncol(data_eset1), refeset=data_eset)

  try(arrayQualityMetrics::arrayQualityMetrics(expressionset=data_eset1,
                                               outdir=paste0(dataname,"_post_norm_QC_report"),
                                               force=TRUE, do.logtransform=FALSE,
                                               intgroup=Biobase::phenoData(data_eset1)@varMetadata$labelDescription,
                                               reporttitle=paste0(dataname, " Post Normalization QC report")))

  if (todf) {
    data_eset_df <- ms_eSet_to_df(data_eset1, columnorder=names(data))
    if (length(attr(data_eset_df,"outdir"))==0 & length(outdir)>0) {
      attr(data_eset_df,"outdir") <- outdir
    }
    ms_filewrite(data_eset_df, paste0(dataname, "_data_post_normalization.txt"),
                            outdir=outdir)
    return(data_eset_df)
  } else {
    if (length(attr(data_eset1,"outdir"))==0 & length(outdir)>0) {
      attr(data_eset1,"outdir") <- outdir
    }
    return(data_eset1)
  }
}
