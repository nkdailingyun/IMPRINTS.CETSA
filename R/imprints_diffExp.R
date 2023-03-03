#' imprints_diffExp
#'
#' Function to have an overview of the IMPRINTS-CETSA dataset, regarding to
#' expression level changes only, to select out the differentially expressed proteins
#'
#' @param data dataset after imprints_normalization() function, readings in log2 format
#' @param set a single character to indicate the sample name
#' @param basetemp the character indicating the baseline temperature for
#' expression levels, default value is 37C
#' @param contrast a character to indicate the contrasting treatment conditions
#' @param logFC_threshold the threshold value for log fold changes, default set at 0.3
#' @param adjp_threshold the threshold value for adjusted p values, default set at 0.01
#' @param labelnodes whether to label the proteins with significant differential expression
#' @param xlimit a two numeric element vector to indicate the limit of x-axis
#' @param ylimit a two numeric element vector to indicate the limit of y-axis
#' @param labelgeneid a vector of the gene symbol id to show on the plot

#'
#' @importFrom limma contrasts.fit eBayes lmFit makeContrasts topTable
#' @importFrom ggplot2 ggsave
#' @import ggrepel
#' @export
#' @return a dataframe
#' @examples \dontrun{
#'   MOLM <- imprints_diffExp(MOLM, set="M13", contrast="TNFa-DMSO")
#' }
#'
imprints_diffExp <- function(data, set=NULL, basetemp="37C", contrast=NULL,
                          logFC_threshold=0.3, adjp_threshold=0.01, pfdatabase=FALSE,
                          labelnodes=TRUE, xlimit=NULL, ylimit=NULL,
                          labelgeneid=NULL, returneset=FALSE) {

  dataname <- deparse(substitute(data))
  outdir <- ms_directory(data, dataname)$outdir
  data <- ms_directory(data, dataname)$data

  data_a <- imprints_average(data)

  if (length(grep("description", names(data)))) {
    proteininfo <- unique(data[ ,c("id","description")])
  }
  if (length(grep("countNum", names(data)))) {
    countinfo <- unique(data[ ,stringr::str_which(names(data), "^id$|^sumPSM|^countNum|^sumUniPeps")])
    data <- data[ ,-stringr::str_which(names(data), "^sumPSM|^countNum|^sumUniPeps")]
  }

  if (length(set)>0) {
    data <- data[ ,c(1,2,grep(set,names(data)))]
    data_a <- data_a[ ,c(1,2,grep(set,names(data_a)))]
  }

  if (sum(grepl(basetemp,names(data))) & sum(grepl(basetemp,names(data_a)))) {
    data <- data[ ,c(1,2,grep(basetemp,names(data)))]
    data_a <- data_a[ ,c(1,2,grep(basetemp,names(data_a)))]
  } else {stop("Make sure the basetemp info is embeded in the column names.")}

  cname <- setdiff(names(data), c("id","description","sumUniPeps","sumPSMs","countNum"))
  if (length(unlist(strsplit(cname[1], "_")))==3) {
    temperature <- unlist(lapply(strsplit(cname, "_"),`[`,1))
    replicate <- unlist(lapply(strsplit(cname, "_"),`[`,2))
    treatment <- unlist(lapply(strsplit(cname, "_"),`[`,3))
    pdata1 <- data.frame(temperature=temperature, replicate=replicate, treatment=treatment)
    row.names(pdata1) <- cname
    nread <- nrow(pdata1)
    #print(pdata1)
    data <- merge(data, countinfo)
    data_eset <- ms_to_eSet(data=data, nread=nread, pdata=pdata1)
  } else if (length(unlist(strsplit(cname[1], "_")))==4) {
    set <- unlist(lapply(strsplit(cname, "_"),`[`,1))
    temperature <- unlist(lapply(strsplit(cname, "_"),`[`,2))
    replicate <- unlist(lapply(strsplit(cname, "_"),`[`,3))
    treatment <- unlist(lapply(strsplit(cname, "_"),`[`,4))
    pdata1 <- data.frame(set=set, temperature=temperature, replicate=replicate, treatment=treatment)
    row.names(pdata1) <- cname
    nread <- nrow(pdata1)
    #print(pdata1)
    data <- merge(data, countinfo)
    data_eset <- ms_to_eSet(data=data, nread=nread, pdata=pdata1)
  } else {
    stop("make sure the namings of the columns of the dasaset are correct.")
  }

  if (returneset) {return(data_eset)}

  # ct <- factor(unique(pdata1$treatment))
  ct <- factor(pdata1$treatment)
  # print(ct)
  design <- model.matrix(~0+ct)
  colnames(design) <- levels(ct)
  # print(paste0("The design as follow: "))
  # print(design)

  if (length(contrast)) {
    contrast.matrix <- limma::makeContrasts(contrasts=contrast, levels=design)
    # print(contrast.matrix)
  } else {
    stop("pls specify a contrast expression, such as 'TNFa-DMSO'.")
  }
  fit <- limma::lmFit(data_eset, design)
  fit1 <- limma::contrasts.fit(fit, contrast.matrix)
  fit2 <- limma::eBayes(fit1)#, proportion=0.2, trend=T)
  fittoptable <- NULL
  for (i in 1:length(contrast)) {
    contrast.name <- colnames(fit2$contrasts)[i]
    top <- limma::topTable(fit2, coef=i, number=Inf, adjust="BH", sort.by="p")
    top <- tibble::rownames_to_column(top, "id")
    top$contrast <- contrast.name
    top$category <- ifelse(abs(top$logFC)>=logFC_threshold & top$adj.P.Val<=adjp_threshold, "C", "N")
    top$category <- factor(top$category, levels=c("C","N"))
    message(paste0("The categories of expression level change in ", contrast.name, " are as follows: "))
    print(table(top$category))
    ms_filewrite(top, paste0(dataname, "_",contrast.name,"_eBays.txt"), outdir=outdir)

    top <- top %>% dplyr::rowwise() %>%
      dplyr::mutate(gene=getGeneName(description,pfdatabase),
                    log10p=-log10(adj.P.Val)) %>%
      dplyr::ungroup()

    if(!length(xlimit)) {xlimit <- c(-max(abs(top$logFC),na.rm=T)-0.1, max(abs(top$logFC),na.rm=T)+0.1)}
    if(!length(ylimit)) {ylimit <- c(min(top$log10p,na.rm=T), max(top$log10p,na.rm=T)+0.1)}
    q <- ggpubr::ggscatter(top, x = "logFC", y = "log10p",
                           color = "category", shape=20, alpha=0.3,
                           palette = c("#FC4E07", "gray"),
                           title = paste0("Expression Changes of ",contrast.name),
                           xlab = "fold change of 37C expression level [log2]",
                           ylab = "adjusted p values [-log10]",
                           xlim=xlimit, ylim=ylimit)
    if (labelnodes) {
      if (length(labelgeneid)) {
        q <- q + ggrepel::geom_text_repel(data=subset(top, gene %in% labelgeneid),
                                          aes(label=gene))
      } else {
        q <- q + ggrepel::geom_text_repel(data=subset(top, category=="C"),
                                          aes(label=gene))
      }
    }
    q <- q + geom_hline(yintercept=-log10(adjp_threshold), linetype="dashed") +
      theme(text=element_text(size=12), plot.title=element_text(hjust=0.5, size=rel(1.2)), aspect.ratio=1)
    ggplot2::ggsave(filename=paste0(outdir,"/",format(Sys.time(), "%y%m%d_%H%M"),
                                    "_Protein_abundance_Changes_in_", contrast.name, ".pdf"), q, width=8.27, height=11.69)

    # to make scatter plot
    name1 <- unlist(strsplit(contrast.name, "-"))[1]
    name2 <- unlist(strsplit(contrast.name, "-"))[2]
    col1 <- grep(paste0("_",name1), names(data_a))
    col2 <- grep(paste0("_",name2), names(data_a))
    data_a1 <- data_a[ ,c(1,2,col1,col2)]
    names(data_a1)[c(3,4)] <- c(name1,name2)
    data_a1 <- merge(data_a1, top[ ,c("id","category","gene")], by="id")
    qs <- ggpubr::ggscatter(data_a1, x = name1, y = name2,
                            color = "category", shape=20, alpha=0.3,
                            palette = c("#FC4E07", "gray"),
                            title = "Protein abundance scatter plot",
                            xlab = paste0("Protein level in ", name1),
                            ylab = paste0("Protein level in ", name2))

    axismin <- min(min(data_a1[,3],na.rm=T), min(data_a1[,4],na.rm=T))-0.1
    axismax <- max(max(data_a1[,3],na.rm=T), max(data_a1[,4],na.rm=T))+0.1
    qs <- qs + coord_cartesian(xlim=c(axismin,axismax), ylim=c(axismin,axismax))
    if (labelnodes) {
      if (length(labelgeneid)) {
        qs <- qs + ggrepel::geom_text_repel(data=subset(data_a1, gene %in% labelgeneid),
                                          aes(label=gene))
      } else {
        qs <- qs + ggrepel::geom_text_repel(data=subset(data_a1, category=="C"),
                                          aes(label=gene))
      }
    }
    qs <- qs + geom_abline(slope=1,intercept=0, linetype="dashed", color="blue") +
      theme(text=element_text(size=12), plot.title=element_text(hjust=0.5, size=rel(1.2)), aspect.ratio=1)
    ggplot2::ggsave(filename=paste0(outdir,"/",format(Sys.time(), "%y%m%d_%H%M"),
                                    "_Protein_abundance_in_", contrast.name, ".pdf"), qs, width=8.27, height=11.69)

    fittoptable <- rbind(fittoptable, top)
  }

  if (length(attr(fittoptable,"outdir"))==0 & length(outdir)>0) {
    attr(fittoptable,"outdir") <- outdir
  }
  return(fittoptable)
}
