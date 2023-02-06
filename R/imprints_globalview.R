#' imprints_globalview
#'
#' Function to have a global overview of the IMPRINTS-CETSA dataset, regarding to
#' protein abundance level changes and stability shifts
#'
#' @param data dataset after imprints_caldiff() function
#' @param basetemp the character indicating the baseline temperature for
#' abundance/expression levels, default value is 37C
#' @param qc whether to apply the reproducibility quality control on input data,
#' default set to TRUE
#' @param cvthreshold the CV threshold value for subsetting reproducible
#' measurements, default value is 0.1
#' @param corrthreshold the Correlation threshold value for subsetting
#' reproducible measurements, default value is 0.5
#' @param abundancechange_nMAD the number of MADs to set the significance
#' cutoff on protein abundance level change, default value is 2.5
#' @param stabilitychange_nMAD the number of MADs to set the significance cutoff
#' on protein thermal stability shift, default value is 2.5
#' @param labelnodes whether to text-label the selected nodes, default set to FALSE
#' @param labelcategory the categories of nodes to label, default value is c("CC","NC","CN")
#' @param labelgeneid a vector of the gene symbol id to show on the plot, exclusive from labelcategory
#'
#' @import dplyr Biobase
#' @import ggpubr
#' @import ggrepel
#' @export
#' @return a dataframe
#' @examples \dontrun{
#'   MOLM <- imprints_globalview(MOLM)
#' }
#'
imprints_globalview <- function(data, basetemp="37C", qc=TRUE, cvthreshold=0.1, corrthreshold=0.5,
                             abundancechange_cutoff=NULL, stabilitychange_cutoff=NULL,
                             abundancechange_nMAD=2.5, stabilitychange_nMAD=2.5,
                             labelnodes=FALSE, labelcategory=c("CC","NC","CN"), labelgeneid=NULL,
                             xrange=NULL, yrange=NULL) {

  dataname <- deparse(substitute(data))
  outdir <- ms_directory(data, dataname)$outdir
  data <- ms_directory(data, dataname)$data

  if (length(grep("description", names(data)))) {
    proteininfo <- unique(data[ ,c("id","description")])
    data$description <- NULL
  }
  if (length(grep("countNum", names(data)))) {
    countinfo <- unique(data[ ,c("id","sumUniPeps","sumPSMs","countNum")])
    data <- data[ ,!(names(data) %in% c("sumUniPeps","sumPSMs","countNum"))]
  }

  refcol <- which(apply(data[,-1],2,sum,na.rm=T)==0)
  datal <- gather(data[ ,-(refcol+1)], condition, reading, -id)
  if (length(unlist(strsplit(datal$condition[1], "_")))==4) {
    datal1 <- tidyr::separate(datal, condition, into=c("set","temperature","replicate","treatment"), sep="_")
    datal <- datal1 %>% group_by(id, set, treatment, temperature) %>%
      summarize(mreading=mean(reading,na.rm=T), sdreading=sd(reading, na.rm=T))
    datal_copy <- datal
    datal_copy <- merge(proteininfo, datal_copy)
    names(datal_copy)[c(6:7)] <- c("arithmetic mean","standard deviation")
    ms_filewrite(datal_copy, paste0("summaried_readings_in_", dataname, ".txt"), outdir=outdir)
    if (qc) {
      reproducible1 <- datal %>% group_by(id,set,treatment) %>% summarise(cvreading=sqrt(exp(mean(sdreading,na.rm=T)^2)-1)) %>%
        filter(cvreading<=cvthreshold)
      reproducible2 <- tidyr::spread(datal1, replicate, reading)
      reproducible2 <- plyr::ddply(reproducible2, "id", function(data) {
        a <- cor(data[ ,-c(1:4)], use="complete.obs")
        data.frame(corr=mean(a[lower.tri(a)]))
      })
      reproducible2 <- subset(reproducible2, corr>=corrthreshold)
      datal <- subset(datal, id %in% unique(c(reproducible1$id,reproducible2$id)))
    }
    datal_Cchange <- datal %>% group_by(id,set,treatment) %>% summarize(change=max(mreading,na.rm=T)-min(mreading,na.rm=T))
    datal_Echange <- datal[grep(basetemp, datal$temperature), c(1,2,3,5)]
    datal_change <- na.omit(merge(datal_Echange, datal_Cchange))
    names(datal_change)[c(4,5)] <- c("abundancechange", "stabilitychange")
  } else if (length(unlist(strsplit(datal$condition[1], "_")))==3) {
    datal1 <- tidyr::separate(datal, condition, into=c("temperature","replicate","treatment"), sep="_")
    datal <- datal1 %>% group_by(id, treatment, temperature) %>%
      summarize(mreading=mean(reading,na.rm=T), sdreading=sd(reading, na.rm=T))
    datal_copy <- datal
    datal_copy <- merge(proteininfo, datal_copy)
    names(datal_copy)[c(5:6)] <- c("arithmetic mean","standard deviation")
    ms_filewrite(datal_copy, paste0("summaried_readings_in_", dataname, ".txt"), outdir=outdir)
    if (qc) {
      reproducible1 <- datal %>% group_by(id,treatment) %>% summarise(cvreading=sqrt(exp(mean(sdreading,na.rm=T)^2)-1)) %>%
        filter(cvreading<=cvthreshold)
      reproducible2 <- tidyr::spread(datal1, replicate, reading)
      reproducible2 <- plyr::ddply(reproducible2, "id", function(data) {
        a <- cor(data[ ,-c(1:3)], use="complete.obs")
        data.frame(corr=mean(a[lower.tri(a)]))
      })
      reproducible2 <- subset(reproducible2, corr>=corrthreshold)
      datal <- subset(datal, id %in% unique(c(reproducible1$id,reproducible2$id)))
    }
    datal_Cchange <- datal %>% group_by(id,treatment) %>% summarize(change=max(mreading,na.rm=T)-min(mreading,na.rm=T))
    datal_Echange <- datal[grep(basetemp, datal$temperature), c(1,2,4)]
    datal_change <- na.omit(merge(datal_Echange, datal_Cchange))
    names(datal_change)[c(3,4)] <- c("abundancechange", "stabilitychange")
  }

  if (length(abundancechange_cutoff)==1) {
    abundancecutoff <- abundancechange_cutoff
  } else {
    abundancecutoff <- median(datal_change$abundancechange,na.rm=T)+abundancechange_nMAD*mad(datal_change$abundancechange,na.rm=T)
    print(paste0("Abundance level change cutoff set at ", round(abundancecutoff,3)))
  }
  if (length(stabilitychange_cutoff)==1) {
    stabilitycutoff <- stabilitychange_cutoff
  } else {
    stabilitycutoff <- median(datal_change$stabilitychange,na.rm=T)+stabilitychange_nMAD*mad(datal_change$stabilitychange,na.rm=T)
    print(paste0("Stability level change cutoff set at ", round(stabilitycutoff,3)))
  }

  datal_change <- datal_change %>% rowwise() %>%
    mutate(abundance.hit = ifelse(abs(abundancechange)<abundancecutoff, F, T)) %>%
    #mutate(basechangedir = ifelse(mreading < 0, "-", "+")) %>%
    mutate(stability.hit = ifelse(stabilitychange<stabilitycutoff, F, T)) %>%
    mutate(category=paste0(abundance.hit, stability.hit))
  datal_change$category <- gsub("FALSE","N",datal_change$category)
  datal_change$category <- gsub("TRUE","C",datal_change$category)

  datal_change <- proteininfo %>% rowwise() %>% mutate(gene=getGeneName(description, pfdatabase)) %>%
    inner_join(datal_change) %>% arrange(category)

  print(paste0("The category of abundance level change and thermal stability shift are as follows: "))
  print(table(datal_change$category))

  if(!length(xrange)) {
    xrange <- c(-max(abs(datal_change$abundancechange),na.rm=T)-0.1, max(abs(datal_change$abundancechange),na.rm=T)+0.1)
  }
  if(!length(yrange)) {
    yrange <- c(0, max(abs(datal_change$stabilitychange),na.rm=T)+0.1)
  }

  if (length(grep("set",names(datal_change)))) {
    q <- ggpubr::ggscatter(datal_change, x = "abundancechange", y = "stabilitychange",
                           color = "category", shape=20, alpha=0.9, facet.by=c("set","treatment"),
                           palette = c("#FC4E07", "#00AFBB", "#E7B800", "gray"),
                           xlab = "37C protein abundance level change",
                           ylab = "Delta fold change across temperature",
                           xlim = xrange, ylim = yrange)
  } else {
    q <- ggpubr::ggscatter(datal_change, x = "abundancechange", y = "stabilitychange",
                           color = "category", shape=20, alpha=0.9, facet.by="treatment",
                           palette = c("#FC4E07", "#00AFBB", "#E7B800", "gray"),
                           xlab = "37C protein abundance level change",
                           ylab = "Delta fold change across temperature",
                           xlim = xrange, ylim = yrange)
  }

  if (labelnodes) {
    if (length(labelgeneid)) {
      q <- q + ggrepel::geom_text_repel(data=subset(datal_change, gene %in% labelgeneid),
                                        aes(label=gene), arrow=arrow(length = unit(0.02, "npc")),
                                        max.overlaps=50) + geom_point(data=subset(datal_change, gene %in% labelgeneid),color="black")
    } else {
      q <- q + ggrepel::geom_text_repel(data=subset(datal_change, category %in% labelcategory),
                                        aes(label=gene), arrow=arrow(length = unit(0.02, "npc")),
                                        max.overlaps=50)
    }
  }

  q <- q + geom_hline(yintercept=stabilitycutoff, linetype="dashed", color="black") +
    geom_vline(xintercept=-abundancecutoff, linetype="dashed", color="black") +
    geom_vline(xintercept=abundancecutoff, linetype="dashed", color="black") +
    theme(text=element_text(size=12), plot.title=element_text(hjust=0.5, size=rel(1.2)), aspect.ratio=1)

  ggsave(file=paste0(outdir,"/",format(Sys.time(), "%y%m%d_%H%M"), "_GlobaView_of_Changes_in_", dataname, ".pdf"), q, width=11.69, height=8.27)
  ms_filewrite(datal_change, paste0("Changes_in_", dataname, ".txt"), outdir=outdir)

  return(datal_change)
}
