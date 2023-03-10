% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/imprints_score_abundance.R
\name{imprints_score_abundance}
\alias{imprints_score_abundance}
\title{imprints_score_abundance}
\usage{
imprints_score_abundance(
  data,
  set = NULL,
  contrast = NULL,
  basetemp = "37C",
  pfdatabase = FALSE,
  logFC_threshold = 0.2,
  adjp_threshold = 0.01,
  labelnodes = FALSE,
  labelgeneid = NULL,
  returnsplitlist = FALSE
)
}
\arguments{
\item{data}{dataset after imprints_normalization() function, readings in log2 format}

\item{set}{a single character to indicate the sample name}

\item{contrast}{a character to indicate the contrasting treatment conditions}

\item{basetemp}{character indicating the baseline temperature for expression levels, default value is 37C}

\item{logFC_threshold}{the threshold value for log fold changes, default set at 0.3}

\item{adjp_threshold}{the threshold value for adjusted p values, default set at 0.01}

\item{labelnodes}{whether to label the proteins with significant differential expression}

\item{labelgeneid}{a vector of the gene symbol id to show on the plot}
}
\value{
a dataframe
}
\description{
Function to calculate the relative abundance change at each temperature,
based on the normalized abundance data but not the logratio data
}
\examples{
\dontrun{
  MOLM <- imprints_score_abundance(MOLM, set="M13", contrast="TNFa-DMSO")
}

}
