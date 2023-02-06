#' ms_innerbarplotting
#'
#' Internal function to make barplots, with new features from Marc
#' @param data dataset to plot in barplot format
#' @param withset whether the data contains more than one set of data
#' @param witherrorbar a logical to print or not the error bar on the plot
#' @param colorpanel a vector of customizable color scheme provided by the user, default set
#' c("gray","blue","orange")
#' @param usegradient whether the barplot should be draw in color gradient format
#' @param colorgradient the color scheme of gradient applied, default value c("#4575B4","ivory", "#D73027")
#' @param linegraph whether to plot the graph in a line graph format, default set to FALSE
#' @param log2scale whether the yscales should be in log2 scale, default set to TRUE
#' @param ratio aspect ratio of the plot, default set to 0.6
#' @param layout a vector indicating the panel layout for multi-panel plots
#' per page, default value is c(2,3) for set containing data, otherwise c(4,3)
#' @keywords internal

#' @return a command
#' @examples \dontrun{
#' }
#'

ms_innerbarplotting <- function(data, withset=FALSE,
                        witherrorbar=TRUE, colorpanel=c("gray","blue","orange"),
                        usegradient=FALSE, colorgradient=c("#4575B4","ivory", "#D73027"),
                        linegraph=FALSE, log2scale=TRUE, ratio=0.6, layout=NULL) {

  # print(max(data$reading))
  # print(min(data$reading))
  if (withset) {
    data <- droplevels(data)
    data_list <- split(data, data$set)
    q_list <- list()
    n_loop <- 1
    for (j in names(data_list)) {
      # print(data_list[[j]])
      if (nrow(data_list[[j]])>0) {
        d2 <- data_list[[j]]
        if (!log2scale) {
          minreading=0.5
          maxreading=2
          legendscale = c(min(max(min(d2$mean, na.rm=T)-0.5, 0), minreading), max(max(d2$mean, na.rm=T)+0.5, maxreading))
        } else {
          minreading=-0.5
          maxreading=0.5
          legendscale = c(min(min(d2$mean, na.rm=T)-0.1, minreading), max(max(d2$mean, na.rm=T)+0.1, maxreading))
        }
        # print(legendscale)
        q <- ggplot(d2, aes(x=condition, y=mean, fill=treatment)) +
          geom_bar(stat="identity") + coord_cartesian(ylim=legendscale) +
          scale_fill_manual(drop=FALSE, values=colorpanel)
        if (witherrorbar) { q <- q + geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2, # Width of the error bars
                                                   position=position_dodge(.9)) }
        if (log2scale) {
          q <- q + ylab("fold change(log2)") + ggtitle(paste(j, as.character(unique(d2$id)), sep="\n"))
        } else {
          q <- q + ylab("fold change") + ggtitle(paste(j, as.character(unique(d2$id)), sep="\n"))
        }

        q <- q + theme_classic() +
          # labs(subtitle = subt$category[n_loop]) +
          theme(
            text = element_text(size=10),
            strip.text.x = element_text(size = 5),
            plot.title = element_text(hjust=0.5, size = rel(0.8)),
            legend.background=element_rect(fill=NULL),
            legend.key.height=unit(0.5, "cm"),
            legend.key.width=unit(0.15, "cm"),
            legend.title = element_text(face = "bold"),
            legend.text = element_text(size = rel(0.7)),
            legend.justification="center",
            #legend.position="none", #c(0.2,0.8),
            panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            strip.background=element_blank(),
            #strip.text = element_blank(),
            axis.line.x = element_line(),
            axis.line.y = element_line(),
            axis.text.x = element_text(angle = 45, hjust = 1, size = rel(0.7)),
            aspect.ratio = 0.6
          )
        q_list[[j]] <- q
      } else {
        q <- ggplot()
        q_list[[j]] <- q
      }
    }
    q_list <- gridExtra::grid.arrange(grobs=q_list, ncol=1)
    return(q_list)
  } else {
    if (!log2scale) {
      minreading=0.5
      maxreading=2
      legendscale = c(min(max(min(data$mean, na.rm=T)-0.5, 0), minreading), max(max(data$mean, na.rm=T)+0.5, maxreading))
    } else {
      minreading=-0.5
      maxreading=0.5
      legendscale = c(min(min(data$mean, na.rm=T)-0.1, minreading), max(max(data$mean, na.rm=T)+0.1, maxreading))
    }
    if (linegraph) {
      q <- ggplot(data, aes(x=treatment, y=mean, group=temperature, color=temperature)) +
        geom_line() + geom_point() + coord_cartesian(ylim=legendscale) +
        scale_color_manual(drop=FALSE, values=colorpanel)
    } else if (!usegradient) {
      q <- ggplot(data, aes(x=condition, y=mean, fill=treatment)) +
        geom_bar(stat="identity") + coord_cartesian(ylim=legendscale) +
        scale_fill_manual(drop=FALSE, values=colorpanel)
    } else {
      q <- ggplot(data, aes(x=condition, y=mean, fill=mean)) +
        geom_bar(stat="identity") + coord_cartesian(ylim=legendscale) +
        scale_fill_gradient2(limits=legendscale, low=colorgradient[1], mid=colorgradient[2], high=colorgradient[3],
                             midpoint=0, na.value="gray90", guide=guide_colorbar(""))
    }

    if (witherrorbar) {
      if (linegraph) {
        q <- q + geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.1) # Width of the error bars
      } else {
        q <- q + geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2, # Width of the error bars
                               position=position_dodge(.9))
      }
    }
    if (log2scale) {
      q <- q + ylab("fold change(log2)") + ggtitle(as.character(unique(data$id)))
    } else {
      q <- q + ylab("fold change") + ggtitle(as.character(unique(data$id)))
    }

    q <- q + theme_classic() +
      theme(
        text = element_text(size=10),
        strip.text.x = element_text(size = 5),
        plot.title = element_text(hjust=0.5, size = rel(0.8)),
        legend.background=element_rect(fill=NULL),
        legend.key.height=unit(0.5, "cm"),
        legend.key.width=unit(0.15, "cm"),
        legend.title = element_text(face = "bold"),
        legend.text = element_text(size = rel(0.7)),
        legend.justification="center",
        #legend.position="none", #c(0.2,0.8),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        strip.background=element_blank(),
        #strip.text = element_blank(),
        axis.line.x = element_line(),
        axis.line.y = element_line(),
        axis.text.x = element_text(angle = 45, hjust = 1, size = rel(0.7)),
        aspect.ratio = ratio
      )
    return(q)
  }
}
