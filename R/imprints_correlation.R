#' imprints_correlation
#'
#' Function to calculate the Pearson correlation of two treatment conditions
#'
#' @param data dataset of average IMPRINTS profile, that is after calculating the relative protein
#' abundance differences and deriving the average, i.e., after imprints_average(), make sure
#' the columns with readings are named in the format like "37C_Treatment"
#' @param setvector if necessary, a two-element character vector with distinguishable keyword
#' from the column names of the input dataset
#' @param treatmentvector a two-element character vector with distinguishable keyword
#' from the column names of the input dataset
#' @param bygene whether to match the protein accessions by gene, default set to TRUE,
#' otherwise would match be protein accession which could be affected by protein isoforms,
#' when there is isoform difference, the first dataset would be used as the reference name
#'
#'
#' @import dplyr
#' @export
#' @return a dataframe
#' @examples \dontrun{
#'   MOLM <- imprints_correlation(MOLM, setvector=NULL, treatmentvector=c("T1","T2"))
#'   MOLM <- imprints_correlation(MOLM, setvector=c("M13","M16"), treatmentvector=c("T1","T2"))
#' }
#'
#'

imprints_correlation <- function(data=NULL, setvector=NULL, treatmentvector=NULL, bygene=TRUE) {

  if (!length(setvector) %in% c(0,2)) {stop("Only correlation of two treatment were allowed,
                                       the number of set could be either 0 or 2")}
  if (length(treatmentvector) != 2) {stop("Only correlation of two treatment were allowed,
                                    pls specify a two-element vector with distinguishable keyword")}

  if (inherits(data,"data.frame") | class(data)=="data.frame") {
    dataname <- deparse(substitute(data))
    outdir <- ms_directory(data, dataname)$outdir
    data <- ms_directory(data, dataname)$data

    subset1 <- grep(paste0("_",treatmentvector[1]),names(data))
    subset2 <- grep(paste0("_",treatmentvector[2]),names(data))
    if (length(setvector)) {
      subset1_set <- grep(setvector[1],names(data))
      subset2_set <- grep(setvector[2],names(data))
      subset1 <- intersect(subset1,subset1_set)
      subset2 <- intersect(subset2,subset2_set)
    }
    if (length(subset1)>0 & length(subset2)>0 & length(subset1)==length(subset2)) {
      data1 = data[ ,c(1,2,subset1)]
      data2 = data[ ,c(1,2,subset2)]
    } else {
      stop("pls specify a two-element vector with distinguishable keyword from input dataset")
    }
  } else if (class(data)=="list" & length(data)==2) {
    dataname <- deparse(substitute(data[[1]]))
    outdir <- ms_directory(data[[1]], dataname)$outdir

    subset1 <- grep(treatmentvector[1],names(data[[1]]))
    subset2 <- grep(treatmentvector[2],names(data[[2]]))
    if (length(setvector)) {
      subset1_set <- grep(setvector[1],names(data[[1]]))
      subset2_set <- grep(setvector[2],names(data[[2]]))
      subset1 <- intersect(subset1,subset1_set)
      subset2 <- intersect(subset2,subset2_set)
    }
    if (length(subset1)>0 & length(subset2)>0 & length(subset1)==length(subset2)) {
      data1 = data[[1]][ ,c(1,2,subset1)]
      data2 = data[[2]][ ,c(1,2,subset2)]
    } else {
      stop("pls specify a two-element vector with distinguishable keyword from input dataset")
    }
  }

  if (length(grep("countNum", names(data1)))) {
    countinfo1 <- unique(data1[ ,c("id","sumUniPeps","sumPSMs","countNum")])
    data1 <- data1[ ,!(names(data1) %in% c("sumUniPeps","sumPSMs","countNum"))]
  }
  if (length(grep("countNum", names(data2)))) {
    countinfo2 <- unique(data2[ ,c("id","sumUniPeps","sumPSMs","countNum")])
    data2 <- data2[ ,!(names(data2) %in% c("sumUniPeps","sumPSMs","countNum"))]
  }

  if (bygene) {
    # merge by the gene symbol
    data1 <- data1 %>% rowwise() %>% mutate(gene = getGeneName(description))
    if (anyNA(data1$gene)) {
      data1 <- data1[-which(is.na(data1$gene)), ]
    }
    data2 <- data2 %>% rowwise() %>% mutate(gene = getGeneName(description))
    if (anyNA(data2$gene)) {
      data2 <- data2[-which(is.na(data2$gene)), ]
    }
    proteininfo1 <- data1[ ,c("id","description","gene")]
    proteininfo2 <- data2[ ,c("id","description","gene")]
    # print(head(data1))
    # print(head(data2))
    data <- merge(data1[ ,-c(1:2)], data2[ ,-c(1:2)], by="gene")
    # print(head(data))
  } else {
    # merge by the Uniprot ID
    proteininfo1 <- data1[ ,c("id","description")]
    proteininfo2 <- data2[ ,c("id","description")]
    # print(head(data1))
    # print(head(data2))
    data <- merge(data1[ ,-2], data2[ ,-2], by="id")
    # print(head(data))
  }
  # return(data)
  print(paste0("The number of protein entries in common is ", nrow(data)))
  #print(unique(data$gene))
  if (bygene & sum(duplicated(data$gene))) {
    print(paste0(sum(duplicated(data$gene)), " genes were duplicated and removed: "))
    dup_gene <- sort(unique(data$gene[duplicated(data$gene)]))
    print(dup_gene)
    data <- data[-which(data$gene %in% dup_gene), ]
  }
  # return(data)
  if (bygene) {
    datal <- gather(data, treatment, reading, -gene)
  } else {
    datal <- gather(data, treatment, reading, -id)
  }
  if (length(unlist(strsplit(datal$treatment[1], "_")))==3) {
    datal <- separate(datal, treatment, into=c("set","temperature","condition"))
    datal <- unite(datal, "condition", c("set","condition"))
  } else if (length(unlist(strsplit(datal$treatment[1], "_")))==2) {
    datal <- separate(datal, treatment, into=c("temperature","condition"))
  }
  treatment = unique(datal$condition)
  print(paste0("The conditions for correlation calculation are ", treatment[1], " and ", treatment[2]))
  dataw <- spread(datal, condition, reading)
  # print(as.character(eval(parse(text = treatment[1]))))
  # print(head(dataw))
  # return(list(datal=datal, dataw=dataw, treatment=treatment))
  if (bygene) {
    cortable <- group_by(dataw,gene) %>%
      summarize(correlation=cor(eval(parse(text=treatment[1])),eval(parse(text=treatment[2])), use="complete.obs"))
    dataw <- spread(datal, temperature, reading)
    dism <- plyr::ddply(dataw, "gene", function(data) {
      data<-data[order(data$condition), ]
      dm<-as.matrix(dist(data[ ,-c(1:2)]))[1,2]
    })
    names(dism)[2] <- "distance"
    # print(head(cortable))
    # hist(cortable$correlation, breaks=100, xlab="", main="Distribution of Pearson correlations")
    cortable <- merge(proteininfo1,cortable)
    cortable <- merge(cortable,dism)
    return(cortable[order(cortable$correlation,decreasing=T),-1])
  } else {
    cortable <- group_by(dataw,id) %>%
      summarize(correlation=cor(eval(parse(text=treatment[1])),eval(parse(text=treatment[2])), use="complete.obs"))
    dataw <- spread(datal, temperature, reading)
    dism <- plyr::ddply(dataw, "id", function(data) {
      data<-data[order(data$condition), ]
      dm<-as.matrix(dist(data[ ,-c(1:2)]))[1,2]
    })
    names(dism)[2] <- "distance"
    # print(head(cortable))
    # hist(cortable$correlation, breaks=100, xlab="", main="Distribution of Pearson correlations")
    cortable <- merge(proteininfo1,cortable)
    cortable <- merge(cortable,dism)
    cortable <- cortable[order(cortable$correlation,decreasing=T), ]

    if (length(attr(cortable,"outdir"))==0 & length(outdir)>0) {
      attr(cortable,"outdir") <- outdir
    }
    return(cortable)
  }
}
