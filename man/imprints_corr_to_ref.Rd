% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/imprints_corr_to_ref.R
\name{imprints_corr_to_ref}
\alias{imprints_corr_to_ref}
\title{imprints_corr_to_ref}
\usage{
imprints_corr_to_ref(
  data = NULL,
  set = NULL,
  treatment = NULL,
  reference = NULL,
  use_score = c("euclidean", "pearson"),
  score_threshold = 0.9,
  include_neg = FALSE,
  max_na = 0
)
}
\arguments{
\item{data}{dataset of average IMPRINTS profile, that is after calculating the relative protein
abundance differences and deriving the average, i.e., after imprints_average(), make sure
the columns with readings are named in the format like "(Set_)37C_Treatment"}

\item{set}{a single character to indicate the sample name}

\item{treatment}{a single character to indicate the sample name}

\item{reference}{a numeric vector with the profile readings}

\item{use_score}{a single character element that define the method score. Method available : 'euclidean' or 'pearson'}

\item{score_threshold}{a numeric value to indicate the correlation score threshold, default set to 0.9}

\item{include_neg}{whether to include the negatively correlated proteins, default set to FALSE}

\item{max_na}{an integer indicating the maximum number of missing value for one protein, default is 0}
}
\value{
a dataframe with the correlation from the profile and the proteins which passed the threshold
}
\description{
Function to find out the proteins with most similar profile to the reference profile
}
\examples{
\dontrun{
MOLM_M13_T1 <- imprints_corr_to_ref(MOLM, set="M13", treatment="T1",
  reference=c(0.02792811,0.03133724,0.14457743,0.33304728,0.27218208,0.23847792))
}


}
