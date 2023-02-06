#' imprints_score
#'
#' Function to calculate the relative abundance score and thermal stability score
#' Based on the scores, it divides proteins into 9 or 4 categories and shows in a scatter plot
#'
#' @param data dataset after calculating the relative protein abundance differences
#' i.e., imprints_caldiff()
#' @param format choose between 9 and 4, indicating how many categories to segregate,
#' default value is 9, in which proteins are divided into 9 categories: NN, CN+, CN-,
#' NC+, NC-, CC++, CC+-, CC-+, CC--. When switch to 4, proteins are instead divided into
#' 4 categories: NN, CN, NC, CC. The sign of + or - after N or C is determined by the
#' sign of abundance.score and stability.score.mean, respectively
#' @param basetemp the character indicating the baseline temperature for
#' expression levels, default value is 37C
#' @param basecvcheck whether to perform a quality check on the reproducibility of baseline
#' protein expression levels across the replicates, default set to FALSE
#' @param stabilitycvcheck whether to perform a quality check on the reproducibility of protein
#' stability trend across the replicates, default set to FALSE
#' @param cvcutoffthreshold the significance level of threshold used for CV quality control,
#' default value is 2.5, typically can range from 2 to 3
#' @param allzscore when calculating the z score, whether to use all the readings from
#' all the different treatment groups, default set to TRUE
#' @param weightbycv whether to weight the thermal shift with the cv information, default set to TRUE
#' @param fdrthreshold the significance level of global fdr, default value is 0.01
#' @param useMAD whether to use MAD scheme for significance analysis, default set to FALSE
#' @param nMAD the significance level of MAD, default value 2.5
#' @param cutoffvector a vector to specify the symmetrical cutoff values (abundance, stability)
#' @param labelnodes whether to text-label the selected nodes, default set to FALSE
#' @param labelcategory the categories of nodes to label, default value is c("CC","NC","CN")
#' @param labelgeneid a vector of the gene symbol id to show on the plot, exclusive from labelcategory
#'
#' @import dplyr fdrtool
#' @export
#' @return a dataframe
#' @examples \dontrun{
#'   MOLM <- imprints_score(MOLM)
#' }
#'
#'


imprints_score <- function(data, format="9", basetemp="37C", basecvcheck=FALSE, stabilitycvcheck=FALSE,
                           pfdatabase=FALSE, cvcutoffthreshold=2.5, allzscore=TRUE, weightbycv=TRUE,
                           fdrthreshold=0.01, useMAD=FALSE, nMAD=2.5, cutoffvector=NULL, labelnodes=FALSE,
                           labelcategory=c("CC","NC","CN"), labelgeneid=NULL) {

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
    countinfo <- unique(data[ ,c("id","sumUniPeps","sumPSMs","countNum")])
    data <- data[ ,!(names(data) %in% c("sumUniPeps","sumPSMs","countNum"))]
  }

  refcol <- which(apply(data[,-1],2,sum,na.rm=T)==0)
  data <- gather(data[ ,-(refcol+1)], condition, reading, -id)
  if (length(unlist(strsplit(data$condition[1], "_")))==4) {
    data <- tidyr::separate(data, condition, into=c("set","temperature","replicate","treatment"), sep="_")
    data_abd <- data %>% group_by(id, set, treatment) %>%
      filter(temperature==basetemp) %>%
      summarise(abundance.score=mean(reading,na.rm=T),
                abundance.score.cv=sd(reading,na.rm=T)/exp(abs(abundance.score)))
    baseCVcutoff <- mean(data_abd$abundance.score.cv,na.rm=T)+
      cvcutoffthreshold*sd(data_abd$abundance.score.cv,na.rm=T)
    data_thermal <- data %>% left_join(data_abd) %>% rowwise() %>%
      mutate(stability=reading-abundance.score) %>% filter(temperature!=basetemp)
    data_thermal1 <- data_thermal %>% group_by(id,set,treatment,temperature) %>%
      summarise(stability.score=mean(stability,na.rm=T),
                stability.score.cv=sd(stability,na.rm=T)/exp(abs(stability.score)))
    data_thermal1t <- na.omit(data_thermal) %>% group_by(id,set,treatment,temperature) %>%
      summarise(stability.score.t=t.test(stability)$statistic)
    data_thermal1 <- merge(data_thermal1, data_thermal1t, all.x=T)
    if (weightbycv) {
      data_thermal2 <- data_thermal1 %>% group_by(id,set,treatment) %>%
        summarise(stability.score.mean=weighted.mean(stability.score,log(1/stability.score.cv),na.rm=T))
    } else {
      data_thermal2 <- data_thermal1 %>% group_by(id,set,treatment) %>%
        summarise(stability.score.mean=mean(stability.score,na.rm=T))
    }
    data_thermal1.score <- mutate(data_thermal1, temperature=paste("stability.score",temperature,sep="."))
    data_thermal1.score <- tidyr::spread(data_thermal1.score[,-c(6:7)], temperature, stability.score)
    data_thermal1.t <- mutate(data_thermal1, temperature=paste0("stability.score.",temperature,".t"))
    data_thermal1.t <- tidyr::spread(data_thermal1.t[,-c(5:6)], temperature, stability.score.t)
    data_thermal1.cv <- mutate(data_thermal1, temperature=paste0("stability.score.",temperature,".cv"))
    data_thermal1.cv.cutoff <- mean(data_thermal1.cv$stability.score.cv,na.rm=T)+
      cvcutoffthreshold*sd(data_thermal1.cv$stability.score.cv,na.rm=T)
    data_thermal1.cv <- tidyr::spread(data_thermal1.cv[,-c(5,7)], temperature, stability.score.cv)
    data_thermal3 <- merge(data_thermal1.score, data_thermal1.t)
    data_thermal3 <- merge(data_thermal3, data_thermal1.cv)
    data_thermal3 <- merge(data_thermal3, data_thermal2)
    data_score_all <- merge(data_abd, data_thermal3)
    data_score <- data_score_all
    #data1 <- tidyr::unite(data1, condition, set, temperature, replicate, treatment, sep="_")
  } else if (length(unlist(strsplit(data$condition[1], "_")))==3) {
    data <- tidyr::separate(data, condition, into=c("temperature","replicate","treatment"), sep="_")
    data_abd <- data %>% group_by(id, treatment) %>%
      filter(temperature==basetemp) %>%
      summarise(abundance.score=mean(reading,na.rm=T),
                abundance.score.cv=sd(reading,na.rm=T)/exp(abs(abundance.score)))
    baseCVcutoff <- mean(data_abd$abundance.score.cv,na.rm=T)+
      cvcutoffthreshold*sd(data_abd$abundance.score.cv,na.rm=T)
    data_thermal <- data %>% left_join(data_abd) %>% rowwise() %>%
      mutate(stability=reading-abundance.score) %>% filter(temperature!=basetemp)
    data_thermal1 <- data_thermal %>% group_by(id,treatment,temperature) %>%
      summarise(stability.score=mean(stability,na.rm=T),
                stability.score.cv=sd(stability,na.rm=T)/exp(abs(stability.score)))
    data_thermal1t <- na.omit(data_thermal) %>% group_by(id,treatment,temperature) %>%
      summarise(stability.score.t=t.test(stability)$statistic)
    data_thermal1 <- merge(data_thermal1, data_thermal1t, all.x=T)
    if (weightbycv) {
      data_thermal2 <- data_thermal1 %>% group_by(id,treatment) %>%
        summarise(stability.score.mean=weighted.mean(stability.score,log(1/stability.score.cv),na.rm=T))
      # 1/exp(stability.score.cv) could be another option
    } else {
      data_thermal2 <- data_thermal1 %>% group_by(id,treatment) %>% arrange(temperature) %>%
        summarise(stability.score.mean=mean(stability.score,na.rm=T))
    }
    data_thermal1.score <- mutate(data_thermal1, temperature=paste("stability.score",temperature,sep="."))
    data_thermal1.score <- tidyr::spread(data_thermal1.score[,-c(5:6)], temperature, stability.score)
    data_thermal1.t <- mutate(data_thermal1, temperature=paste0("stability.score.",temperature,".t"))
    data_thermal1.t <- tidyr::spread(data_thermal1.t[,-c(4,5)], temperature, stability.score.t)
    data_thermal1.cv <- mutate(data_thermal1, temperature=paste0("stability.score.",temperature,".cv"))
    data_thermal1.cv.cutoff <- mean(data_thermal1.cv$stability.score.cv,na.rm=T)+
      cvcutoffthreshold*sd(data_thermal1.cv$stability.score.cv,na.rm=T)
    data_thermal1.cv <- tidyr::spread(data_thermal1.cv[,-c(4,6)], temperature, stability.score.cv)
    data_thermal3 <- merge(data_thermal1.score, data_thermal1.t)
    data_thermal3 <- merge(data_thermal3, data_thermal1.cv)
    data_thermal3 <- merge(data_thermal3, data_thermal2)
    data_score_all <- merge(data_abd, data_thermal3)
    data_score <- data_score_all
    #data_score <- tidyr::unite(data_score, condition, temperature, replicate, treatment, sep="_")
  } else {
    stop("make sure the namings of the columns of the dasaset are correct.")
  }

  if (basecvcheck) {
    data_abd$baseoutlier <- data_abd$abundance.score.cv>baseCVcutoff
    questionid1 <- unique(subset(data_score, abundance.score.cv>baseCVcutoff)$id)
    print(paste0("There are ", length(questionid1), " proteins with questionable baseline variance..."))
    question1 <- subset(data_score_all, id%in%questionid1)
    question1 <- merge(question1, data_abd[ ,c("id","treatment","baseoutlier")])
    question1 <- merge(proteininfo, question1)
    ms_filewrite(question1, paste0(dataname, "_proteins_base_variance.txt"), outdir=outdir)
    # print(head(question1))
    baseokid <- setdiff(data_score$id,questionid1)
    data_score <- subset(data_score, id%in%baseokid)
  }
  if (stabilitycvcheck) {
    data_thermal1.cv$numberofcvoutlier <- apply(data_thermal1.cv[,-c(1:2)], 1, function(x) sum(x>data_thermal1.cv.cutoff))
    questionid2 <- unique(subset(data_thermal1.cv, numberofcvoutlier>(ncol(data_thermal1.cv)-3)/2)$id)
    print(paste0("There are ", length(questionid2), " proteins with questionable stability variance..."))
    question2 <- subset(data_score_all, id%in%questionid2)
    question2 <- merge(question2, data_thermal1.cv[ ,c("id","treatment","numberofcvoutlier")])
    question2 <- merge(proteininfo, question2)
    ms_filewrite(question2, paste0(dataname, "_proteins_stability_variance.txt"), outdir=outdir)
    # print(head(question1))
    stabilityokid <- setdiff(data_score$id,questionid2)
    data_score <- subset(data_score, id%in%stabilityokid)
  }

  # z-transformation of each score
  data_score_z <- NULL

  if (allzscore) {
    data_score_nocontrol <- data_score[0,]
    for (treat in unique(data_score$treatment)) {
      data_score1 <- subset(data_score, treatment==treat)
      #print(treat)
      if (sum(abs(data_score1$abundance.score),na.rm=T)==0) {next}
      else { data_score_nocontrol <- rbind(data_score_nocontrol, data_score1) }
    }
    col.names <- setdiff(names(data_score_nocontrol),c("id","set","treatment"))
    col.names <- col.names[-grep("\\.cv", col.names)]
    for (cn in col.names) {
      if (grepl("\\.t", cn)) {
        next
      } else {
        ratios <- data_score_nocontrol[ ,cn]
        quants <- quantile(ratios,probs=c(0.1587,0.5,0.8413),na.rm=T)
        ratios[is.na(ratios)] <- 0
        pratios <- ratios>=0
        nratios <- ratios<0
        z <- ratios
        z[nratios] <- ratios[nratios]/as.numeric(diff(quants)[1])
        z[pratios] <- ratios[pratios]/as.numeric(diff(quants)[2])
        data_score_nocontrol[ ,paste0(cn,".z")] <- z
        data_score_nocontrol[ ,paste0(cn,".fdr")] <- fdrtool(z,plot=F,verbose=F)$qval
      }
    }
    data_score_z <- data_score_nocontrol
  } else {
    for (treat in unique(data_score$treatment)) {
      data_score1 <- subset(data_score, treatment==treat)
      #print(treat)
      if (sum(abs(data_score1$abundance.score),na.rm=T)==0) {next}
      col.names <- setdiff(names(data_score1),c("id","treatment"))
      col.names <- col.names[-grep("\\.cv", col.names)]
      for (cn in col.names) {
        if (grepl("\\.t", cn)) {
          next
        } else {
          ratios <- data_score1[ ,cn]
          quants <- quantile(ratios,probs=c(0.1587,0.5,0.8413),na.rm=T)
          ratios[is.na(ratios)] <- 0
          pratios <- ratios>=0
          nratios <- ratios<0
          z <- ratios
          z[nratios] <- ratios[nratios]/as.numeric(diff(quants)[1])
          z[pratios] <- ratios[pratios]/as.numeric(diff(quants)[2])
          data_score1[ ,paste0(cn,".z")] <- z
          data_score1[ ,paste0(cn,".fdr")] <- fdrtool(z,plot=F,verbose=F)$qval
        }
      }
      data_score_z <- rbind(data_score_z, data_score1)
    }
  }

  if (useMAD) {
    if (length(cutoffvector)) {
      cat("Use the user specified cutoff...\n")
      abundancecutoff <- cutoffvector[1]
      print(paste0("Abundance level change cutoff set at ", round(abundancecutoff,3)))
      stabilitycutoff <- cutoffvector[2]
      print(paste0("Stability level change cutoff set at ", round(stabilitycutoff,3)))
    } else {
      abundancecutoff <- median(data_score_z$abundance.score,na.rm=T)+nMAD*mad(data_score_z$abundance.score,na.rm=T)
      print(paste0("Abundance level change cutoff set at ", round(abundancecutoff,3)))
      stabilitycutoff <- median(data_score_z$stability.score.mean,na.rm=T)+nMAD*mad(data_score_z$stability.score.mean,na.rm=T)
      print(paste0("Stability level change cutoff set at ", round(stabilitycutoff,3)))
    }
    if (as.character(format)=="9") {
      data_score_z <- data_score_z %>% rowwise() %>% mutate(abundance.hit=ifelse(abs(abundance.score)>abundancecutoff,T,F))
      data_score_z <- data_score_z %>% rowwise() %>% mutate(stability.hit=ifelse(abs(stability.score.mean)>stabilitycutoff,T,F))
      data_score_z <- data_score_z %>% rowwise() %>% mutate(abundance.sign=ifelse(abundance.score>0,T,F))
      data_score_z <- data_score_z %>% rowwise() %>% mutate(stability.sign=ifelse(stability.score.mean>0,T,F))
      data_score_z <- data_score_z %>% mutate(category=paste0(as.character(abundance.hit),as.character(stability.hit),as.character(abundance.sign),
                                                              as.character(stability.sign)))
      data_score_z$category <- gsub("FALSE","N",data_score_z$category)
      data_score_z$category <- gsub("TRUE","C",data_score_z$category)
      data_score_z$category <- gsub("NN[CN]{2}","NN",data_score_z$category)
      data_score_z$category <- gsub("CNC[CN]{1}","CN+",data_score_z$category)
      data_score_z$category <- gsub("CNN[CN]{1}","CN-",data_score_z$category)
      data_score_z$category <- gsub("NC[CN]{1}N","NC-",data_score_z$category)
      data_score_z$category <- gsub("NC[CN]{1}C","NC+",data_score_z$category)
      data_score_z$category <- gsub("CCCC","CC++",data_score_z$category)
      data_score_z$category <- gsub("CCCN","CC+-",data_score_z$category)
      data_score_z$category <- gsub("CCNC","CC-+",data_score_z$category)
      data_score_z$category <- gsub("CCNN","CC--",data_score_z$category)
    } else if (as.character(format)=="4") {
      data_score_z <- data_score_z %>% rowwise() %>% mutate(abundance.hit=ifelse(abs(abundance.score)>abundancecutoff,T,F))
      data_score_z <- data_score_z %>% rowwise() %>% mutate(stability.hit=ifelse(abs(stability.score.mean)>stabilitycutoff,T,F))
      data_score_z <- data_score_z %>% mutate(category=paste0(as.character(abundance.hit),as.character(stability.hit)))
      data_score_z$category <- gsub("FALSE","N",data_score_z$category)
      data_score_z$category <- gsub("TRUE","C",data_score_z$category)
    }
    data_score <- merge(proteininfo, data_score_z)
    print("The number of proteins in each category is as follows:")
    if (length(grep("set",names(data_score)))) {
      print(table(data_score$category, data_score$treatment, data_score$set, useNA="ifany"))
    } else {
      print(table(data_score$category, data_score$treatment, useNA="ifany"))
    }

    ms_filewrite(data_score, paste0(dataname, "_Abundance_Stability_score_MAD.txt"), outdir=outdir)
    data_score_sim <- data_score[ ,c("id","description","treatment","abundance.score","stability.score.mean","category")]
    data_score_sim <- subset(data_score_sim, category!="NN")
    data_score_sim <- data_score_sim[order(data_score_sim$category), ]
    ms_filewrite(data_score_sim, paste0(dataname, "_Abundance_Stability_score_simple_summary.txt"), outdir=outdir)

    data_score <- data_score %>% rowwise() %>% mutate(gene=getGeneName(description, pfdatabase))
    xlimit <- c(-max(abs(data_score$abundance.score))-0.1, max(abs(data_score$abundance.score))+0.1)
    ylimit <- c(-max(abs(data_score$stability.score.mean))-0.1, max(abs(data_score$stability.score.mean))+0.1)
    palette = c(brewer.pal(length(unique(data_score$category))-1, "Spectral"), "gray")
    if (length(grep("set",names(data_score)))) {
      q <- ggpubr::ggscatter(data_score, x = "abundance.score", y = "stability.score.mean",
                             color = "category", shape=20, alpha=0.9, facet.by=c("set","treatment"),
                             palette = palette,
                             xlab = "Protein Abundance score",
                             ylab = "Protein Thermal stability score",
                             xlim = xlimit, ylim = ylimit)
    } else {
      q <- ggpubr::ggscatter(data_score, x = "abundance.score", y = "stability.score.mean",
                             color = "category", shape=20, alpha=0.9, facet.by="treatment",
                             palette = palette,
                             xlab = "Protein Abundance score",
                             ylab = "Protein Thermal stability score",
                             xlim = xlimit, ylim = ylimit)
    }
    q <- q + geom_hline(yintercept=-stabilitycutoff, linetype="dashed", color = "black") +
      geom_hline(yintercept=stabilitycutoff, linetype="dashed", color = "black") +
      geom_vline(xintercept=-abundancecutoff, linetype="dashed", color = "black") +
      geom_vline(xintercept=abundancecutoff, linetype="dashed", color = "black")

    if (labelnodes) {
      if (length(labelgeneid)) {
        q <- q + ggrepel::geom_text_repel(data=subset(data_score, gene %in% labelgeneid),
                                          aes(label=gene), arrow=arrow(length = unit(0.02, "npc")),
                                          max.overlaps=50) + geom_point(data=subset(data_score, gene %in% labelgeneid))
      } else {
        q <- q + ggrepel::geom_text_repel(data=subset(data_score, category %in% labelcategory),
                                          aes(label=gene), arrow=arrow(length = unit(0.02, "npc")),
                                          max.overlaps=50)
      }
    }
    q <- q + theme(text = element_text(size=12), plot.title = element_text(hjust=0.5, size=rel(1.2)), aspect.ratio=1)
    # q <- q + coord_cartesian(xlim=xrange, ylim=yrange)
    ggsave(file=paste0(outdir,"/",format(Sys.time(), "%y%m%d_%H%M_"),dataname,"_Abundance_Stability_score_MAD.pdf"), q, width=11.69, height=8.27)
  } else {
    if (as.character(format)=="9") {
      data_score_z <- data_score_z %>% rowwise() %>% mutate(abundance.hit=ifelse(abundance.score.fdr<fdrthreshold,T,F))
      data_score_z <- data_score_z %>% rowwise() %>% mutate(stability.hit=ifelse(stability.score.mean.fdr<fdrthreshold,T,F))
      data_score_z <- data_score_z %>% rowwise() %>% mutate(abundance.sign=ifelse(abundance.score>0,T,F))
      data_score_z <- data_score_z %>% rowwise() %>% mutate(stability.sign=ifelse(stability.score.mean>0,T,F))
      data_score_z <- data_score_z %>% mutate(category=paste0(as.character(abundance.hit),as.character(stability.hit),as.character(abundance.sign),
                                                              as.character(stability.sign)))
      data_score_z$category <- gsub("FALSE","N",data_score_z$category)
      data_score_z$category <- gsub("TRUE","C",data_score_z$category)
      data_score_z$category <- gsub("NN[CN]{2}","NN",data_score_z$category)
      data_score_z$category <- gsub("CNC[CN]{1}","CN+",data_score_z$category)
      data_score_z$category <- gsub("CNN[CN]{1}","CN-",data_score_z$category)
      data_score_z$category <- gsub("NC[CN]{1}N","NC-",data_score_z$category)
      data_score_z$category <- gsub("NC[CN]{1}C","NC+",data_score_z$category)
      data_score_z$category <- gsub("CCCC","CC++",data_score_z$category)
      data_score_z$category <- gsub("CCCN","CC+-",data_score_z$category)
      data_score_z$category <- gsub("CCNC","CC-+",data_score_z$category)
      data_score_z$category <- gsub("CCNN","CC--",data_score_z$category)
    } else if (as.character(format)=="4") {
      data_score_z <- data_score_z %>% rowwise() %>% mutate(abundance.hit=ifelse(abundance.score.fdr<fdrthreshold,T,F))
      data_score_z <- data_score_z %>% rowwise() %>% mutate(stability.hit=ifelse(stability.score.mean.fdr<fdrthreshold,T,F))
      data_score_z <- data_score_z %>% mutate(category=paste0(as.character(abundance.hit),as.character(stability.hit)))
      data_score_z$category <- gsub("FALSE","N",data_score_z$category)
      data_score_z$category <- gsub("TRUE","C",data_score_z$category)
    }
    data_score <- merge(proteininfo, data_score_z)
    print("The number of proteins in each category is as follows:")
    if (length(grep("set",names(data_score)))) {
      print(table(data_score$category, data_score$treatment, data_score$set))
    } else {
      print(table(data_score$category, data_score$treatment))
    }
    data_score <- data_score %>% rowwise() %>% mutate(gene=getGeneName(description, pfdatabase)) %>% ungroup()

    ms_filewrite(data_score, paste0(dataname, "_Abundance_Stability_score_z.txt"), outdir=outdir)
    data_score_sim <- data_score[ ,c("id","description","treatment","abundance.score","stability.score.mean","category")]
    data_score_sim <- subset(data_score_sim, category!="NN")
    data_score_sim <- data_score_sim[order(data_score_sim$category), ]
    ms_filewrite(data_score_sim, paste0(dataname, "_Abundance_Stability_score_simple_summary.txt"), outdir=outdir)

    # first to generate the scoring plot on the t-score values
    xlimit <- c(-max(abs(data_score$abundance.score),na.rm=T)-0.1, max(abs(data_score$abundance.score),na.rm=T)+0.1)
    ylimit <- c(-max(abs(data_score$stability.score.mean),na.rm=T)-0.1, max(abs(data_score$stability.score.mean),na.rm=T)+0.1)
    palette = c(brewer.pal(length(unique(data_score$category))-1, "Spectral"), "gray")
    if (length(grep("set",names(data_score)))) {
      q <- ggpubr::ggscatter(data_score, x = "abundance.score", y = "stability.score.mean",
                             color = "category", shape=20, alpha=0.9, facet.by=c("set","treatment"),
                             palette = palette,
                             xlab = "Protein Abundance score",
                             ylab = "Protein Thermal stability score",
                             xlim = xlimit, ylim = ylimit)
    } else {
      q <- ggpubr::ggscatter(data_score, x = "abundance.score", y = "stability.score.mean",
                             color = "category", shape=20, alpha=0.9, facet.by="treatment",
                             palette = palette,
                             xlab = "Protein Abundance score",
                             ylab = "Protein Thermal stability score",
                             xlim = xlimit, ylim = ylimit)
    }

    if (labelnodes) {
      if (length(labelgeneid)) {
        q <- q + ggrepel::geom_text_repel(data=subset(data_score, gene %in% labelgeneid),
                                          aes(label=gene), arrow=arrow(length = unit(0.02, "npc")),
                                          max.overlaps=50) + geom_point(data=subset(data_score, gene %in% labelgeneid))
      } else {
        q <- q + ggrepel::geom_text_repel(data=subset(data_score, category %in% labelcategory),
                                          aes(label=gene), arrow=arrow(length = unit(0.02, "npc")),
                                          max.overlaps=50)
      }
    }
    q <- q + theme(text = element_text(size=12), plot.title = element_text(hjust=0.5, size=rel(1.2)), aspect.ratio=1)
    ggsave(file=paste0(outdir,"/",format(Sys.time(), "%y%m%d_%H%M_"),dataname, "_Abundance_Stability_score.pdf"), q, width=11.69, height=8.27)

    # next to generate the scoring plot on the z-score values
    xlimit <- c(-max(abs(data_score$abundance.score.z),na.rm=T)-0.5, max(abs(data_score$abundance.score.z),na.rm=T)+0.5)
    ylimit <- c(-max(abs(data_score$stability.score.mean.z),na.rm=T)-0.5, max(abs(data_score$stability.score.mean.z),na.rm=T)+0.5)
    if (length(grep("set",names(data_score)))) {
      q <- ggpubr::ggscatter(data_score, x = "abundance.score.z", y = "stability.score.mean.z",
                             color = "category", shape=20, alpha=0.9, facet.by=c("set","treatment"),
                             palette = palette,
                             xlab = "Protein Abundance z-score",
                             ylab = "Protein Thermal stability z-score",
                             xlim = xlimit, ylim = ylimit)
    } else {
      q <- ggpubr::ggscatter(data_score, x = "abundance.score.z", y = "stability.score.mean.z",
                             color = "category", shape=20, alpha=0.9, facet.by="treatment",
                             palette = palette,
                             xlab = "Protein Abundance z-score",
                             ylab = "Protein Thermal stability z-score",
                             xlim = xlimit, ylim = ylimit)
    }
    if (labelnodes) {
      if (length(labelgeneid)) {
        q <- q + ggrepel::geom_text_repel(data=subset(data_score, gene %in% labelgeneid),
                                          aes(label=gene), arrow=arrow(length = unit(0.02, "npc")),
                                          max.overlaps=50) + geom_point(data=subset(data_score, gene %in% labelgeneid))
      } else {
        q <- q + ggrepel::geom_text_repel(data=subset(data_score, category %in% labelcategory),
                                          aes(label=gene), arrow=arrow(length = unit(0.02, "npc")),
                                          max.overlaps=50)
      }
    }
    q <- q + theme(text = element_text(size=12), plot.title = element_text(hjust=0.5, size=rel(1.2)), aspect.ratio=1)
    ggsave(file=paste0(outdir,"/",format(Sys.time(), "%y%m%d_%H%M_"),dataname,"_Abundance_Stability_score_z.pdf"), q, width=11.69, height=8.27)
  }

  if (length(attr(data_score,"outdir"))==0 & length(outdir)>0) {
    attr(data_score,"outdir") <- outdir
  }
  return(data_score)
}
