#' imprints_score_abundance
#'
#' Function to calculate the relative abundance change at each temperature,
#' based on the normalized abundance data but not the logratio data
#'
#' @param data dataset after imprints_normalization() function, readings in log2 format
#' @param set a single character to indicate the sample name
#' @param contrast a character to indicate the contrasting treatment conditions
#' @param basetemp character indicating the baseline temperature for expression levels, default value is 37C
#' @param logFC_threshold the threshold value for log fold changes, default set at 0.3
#' @param adjp_threshold the threshold value for adjusted p values, default set at 0.01
#' @param labelnodes whether to label the proteins with significant differential expression
#' @param labelgeneid a vector of the gene symbol id to show on the plot

#'
#' @import dplyr Biobase
#' @import limma
#' @import ggplot2
#' @import ggrepel
#' @export
#' @return a dataframe
#' @examples \dontrun{
#'   MOLM <- imprints_score_abundance(MOLM, set="M13", contrast="TNFa-DMSO")
#' }
#'
imprints_score_abundance <- function(data, set=NULL, contrast=NULL, basetemp="37C",
                                     pfdatabase=FALSE, logFC_threshold=0.2, adjp_threshold=0.01,
                                     labelnodes=FALSE, labelgeneid=NULL, returnsplitlist=FALSE) {

  dataname <- deparse(substitute(data))
  outdir <- ms_directory(data, dataname)$outdir
  data <- ms_directory(data, dataname)$data

  data_a <- imprints_average(data)

  if (length(grep("description", names(data)))) {
    proteininfo <- unique(data[ ,c("id","description")])
  }
  if (length(grep("countNum", names(data)))) {
    countinfo <- unique(data[ ,c("id","sumUniPeps","sumPSMs","countNum")])
    data <- data[ ,!(names(data) %in% c("sumUniPeps","sumPSMs","countNum"))]
  }

  if (length(set)>0) {
    data <- data[ ,c(1,2,grep(set,names(data)))]
    data_a <- data_a[ ,c(1,2,grep(set,names(data_a)))]
  }

  if (!sum(grepl(basetemp,names(data)))) {
    stop("Make sure the basetemp info is embeded in the column names.")
  }

  cname <- setdiff(names(data), c("id","description","sumUniPeps","sumPSMs","countNum"))
  if (length(unlist(strsplit(cname[1], "_")))==3) {
    temps <- unique(unlist(lapply(strsplit(cname, "_"),`[`,1)))
  } else if (length(unlist(strsplit(cname[1], "_")))==4) {
    temps <- unique(unlist(lapply(strsplit(cname, "_"),`[`,2)))
  }

  fittoptable_total <- NULL
  for (temp in temps) {
    data1 <- data[ ,c(1,2,grep(temp,names(data)))]
    # print(head(data1))
    cname <- setdiff(names(data1), c("id","description","sumUniPeps","sumPSMs","countNum"))
    if (length(unlist(strsplit(cname[1], "_")))==3) {
      temperature <- unlist(lapply(strsplit(cname, "_"),`[`,1))
      replicate <- unlist(lapply(strsplit(cname, "_"),`[`,2))
      treatment <- unlist(lapply(strsplit(cname, "_"),`[`,3))
      pdata1 <- data.frame(temperature=temperature, replicate=replicate, treatment=treatment)
      pdata1 <- mutate(pdata1, trte=paste(treatment,temperature,sep="."))
      row.names(pdata1) <- cname
      nread <- nrow(pdata1)
      #print(pdata1)
      if (!sum(grepl("countNum",names(data1)))) {
        data1 <- merge(data1, countinfo)
      }
      data_eset <- ms_to_eSet(data=data1, nread=nread, pdata=pdata1)
    } else if (length(unlist(strsplit(cname[1], "_")))==4) {
      set <- unlist(lapply(strsplit(cname, "_"),`[`,1))
      temperature <- unlist(lapply(strsplit(cname, "_"),`[`,2))
      replicate <- unlist(lapply(strsplit(cname, "_"),`[`,3))
      treatment <- unlist(lapply(strsplit(cname, "_"),`[`,4))
      pdata1 <- data.frame(set=set, temperature=temperature, replicate=replicate, treatment=treatment)
      pdata1 <- mutate(pdata1, trte=paste(treatment,temperature,sep="."))
      row.names(pdata1) <- cname
      nread <- nrow(pdata1)
      #print(pdata1)
      if (!sum(grepl("countNum",names(data1)))) {
        data1 <- merge(data1, countinfo)
      }
      data_eset <- ms_to_eSet(data=data1, nread=nread, pdata=pdata1)
    } else {
      stop("make sure the namings of the columns of the dasaset are correct.")
    }
    # if (returneset) {return(data_eset)}

    ct <- as.factor(pdata1$trte)
    # print(ct)
    ct2 <- unique(pdata1$temperature)
    design <- model.matrix(~0+ct)
    # print(design)
    colnames(design) <- levels(ct)
    # print(design)

    if (length(contrast)) {
      cons <- c()
      for (i in 1:length(contrast)) {
        con <- strsplit(contrast[i], split="-")[[1]]
        con <- apply(expand.grid(con,ct2), 1, paste, collapse=".")
        con <- paste(con[seq(1,length(con)-1,by=2)],con[seq(2,length(con),by=2)],sep="-")
        cons <- c(cons,con)
      }
      # print(cons)
      contrast.matrix <- makeContrasts(contrasts=cons, levels=design)
      # print(contrast.matrix)
    } else {
      stop("pls specify a contrast expression, such as 'TNFa-DMSO'.")
    }

    # print(dim(data_eset))
    # print(dim(design))
    fit <- lmFit(data_eset, design)
    fit1 <- contrasts.fit(fit, contrast.matrix)
    fit2 <- eBayes(fit1)#, proportion=0.2, trend=T)
    fittoptable <- NULL
    for (i in 1:length(cons)) {
      contrast.name <- colnames(fit2$contrasts)[i]
      top <- topTable(fit2, coef=i, number=Inf, adjust="BH", sort.by="p")
      top <- tibble::rownames_to_column(top, "id")
      top$contrast <- contrast.name
      top$category <- ifelse(abs(top$logFC)>=logFC_threshold & top$adj.P.Val<=adjp_threshold, "C", "N")
      print(paste0("The category of protein change in ", contrast.name, " are as follows: "))
      print(table(top$category,useNA="ifany"))

      top <- top %>% rowwise() %>% mutate(gene=getGeneName(description, pfdatabase), log10p=-log10(adj.P.Val))

      xlimit <- c(-max(abs(top$logFC),na.rm=T)-0.1, max(abs(top$logFC),na.rm=T)+0.1)
      ylimit <- c(0, max(top$log10p,na.rm=T)+0.5)
      q <- ggpubr::ggscatter(top, x = "logFC", y = "log10p",
                             color = "category", shape=20, alpha=0.2,
                             palette = c("#FC4E07", "gray"),
                             title = paste0("Changes of ",contrast.name),
                             xlab = "protein fold change [log2]",
                             ylab = "adjusted p values [-log10]",
                             xlim=xlimit, ylim=ylimit)
      if (labelnodes) {
        if (length(labelgeneid)) {
          q <- q + ggrepel::geom_text_repel(data=subset(top, gene %in% labelgeneid),
                                            aes(label=gene), arrow=arrow(length = unit(0.02, "npc")),
                                            max.overlaps=50)
        } else {
          q <- q + ggrepel::geom_text_repel(data=subset(top, category=="C"),
                                            aes(label=gene), arrow=arrow(length = unit(0.02, "npc")),
                                            max.overlaps=50)
        }
      }
      q <- q + geom_hline(yintercept=-log10(adjp_threshold), linetype="dashed") +
        theme(text=element_text(size=12), plot.title=element_text(hjust=0.5, size=rel(1.2)), aspect.ratio=1)
      ggsave(file=paste0(outdir,"/",format(Sys.time(), "%y%m%d_%H%M"), "_Changes_in_", contrast.name, ".pdf"), q, width=4, height=6)

      fittoptable <- rbind(fittoptable, top)
    }
    fittoptable_total <- rbind(fittoptable_total, fittoptable)
  }

  fittoptable_total <- fittoptable_total[order(fittoptable_total$id), ]
  fittoptable_total$group <- gsub("\\.[0-9]+C","",fittoptable_total$contrast)
  write.csv(fittoptable_total, paste0(outdir,"/",format(Sys.time(),"%y%m%d_%H%M_"),dataname,"_score_by_abundance.csv"), row.names=F)

  q <- ggpubr::ggdensity(fittoptable_total, x="adj.P.Val", fill="lightgray", add="mean", rug=F) +
    facet_wrap(~contrast, ncol=3)
  ggsave(file=paste0(outdir,"/",format(Sys.time(), "%y%m%d_%H%M"), "_adjusted_p_value", ".pdf"), q, width=8.27, height=11.69)

  fittoptable_total1 <- fittoptable_total[grep("37C",fittoptable_total$contrast), ]
  names(fittoptable_total1) <- gsub("category","levelchange",names(fittoptable_total1))

  fittoptable1 <- group_by(fittoptable_total, id, group) %>% summarise(SigNumber=sum(category=="C"))
  fittoptable1 <- merge(unique(fittoptable_total1[,c("id","description","gene","group","levelchange")]),fittoptable1)
  fittoptable1 <- fittoptable1[ ,c("id","description","gene","group","levelchange","SigNumber")]
  write.csv(fittoptable1, paste0(outdir,"/",format(Sys.time(),"%y%m%d_%H%M_"),dataname,"_significant_number.csv"), row.names=F)

  print("The significance test results in the tested conditions are as follows:")
  print(table(fittoptable1$SigNumber, fittoptable1$group))
  if (returnsplitlist) {
    fittoptable1<-split(fittoptable1, f=list(fittoptable1$group,fittoptable1$SigNumber))
  }
  return(fittoptable1)
}
