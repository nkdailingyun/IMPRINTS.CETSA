% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ms_innerbarplotting.R
\name{ms_innerbarplotting}
\alias{ms_innerbarplotting}
\title{ms_innerbarplotting}
\usage{
ms_innerbarplotting(
  data,
  withset = FALSE,
  witherrorbar = TRUE,
  colorpanel = c("gray", "blue", "orange"),
  usegradient = FALSE,
  colorgradient = c("#4575B4", "ivory", "#D73027"),
  linegraph = FALSE,
  log2scale = TRUE,
  ratio = 0.6,
  layout = NULL
)
}
\arguments{
\item{data}{dataset to plot in barplot format}

\item{withset}{whether the data contains more than one set of data}

\item{witherrorbar}{a logical to print or not the error bar on the plot}

\item{colorpanel}{a vector of customizable color scheme provided by the user, default set
c("gray","blue","orange")}

\item{usegradient}{whether the barplot should be draw in color gradient format}

\item{colorgradient}{the color scheme of gradient applied, default value c("#4575B4","ivory", "#D73027")}

\item{linegraph}{whether to plot the graph in a line graph format, default set to FALSE}

\item{log2scale}{whether the yscales should be in log2 scale, default set to TRUE}

\item{ratio}{aspect ratio of the plot, default set to 0.6}

\item{layout}{a vector indicating the panel layout for multi-panel plots
per page, default value is c(2,3) for set containing data, otherwise c(4,3)}
}
\value{
a command
}
\description{
Internal function to make barplots, with new features from Marc
}
\examples{
\dontrun{
}

}
\keyword{internal}
