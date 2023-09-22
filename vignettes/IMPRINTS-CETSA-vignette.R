## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(IMPRINTS.CETSA)

## ---- eval=FALSE--------------------------------------------------------------
#  Chem <- imprints_rawread(c("Chem_37C_3bio_Proteins.txt",
#                             "Chem_47C_3bio_Proteins.txt",
#                             "Chem_50C_3bio_Proteins.txt",
#                             "Chem_52C_3bio_Proteins.txt",
#                             "Chem_54C_3bio_Proteins.txt",
#                             "Chem_57C_3bio_Proteins.txt"),
#                           treatment=c("B1_G1S","B1_S","B1_PM","B2_G1S","B2_S","B2_PM","B3_G1S","B3_S","B3_PM","Mix"))

## ---- eval=FALSE--------------------------------------------------------------
#  Chem <- ms_conditionrename(Chem, incondition = c("37C.1","47C.2","50C.3","52C.4","54C.5","57C.6"),
#                             outcondition=c("37C","47C","50C","52C","54C","57C"))
#  # to remove the "Mix" channel sample
#  Chem <- Chem[,-grep("^Mix$", colnames(Chem))]

## ---- eval=FALSE--------------------------------------------------------------
#  Chem_c <- ms_clean(Chem_c)
#  # to only keep the proteins from the right species, in most cases, this is no more necessary,
#  # because this have been achieved in the ms_clean() function
#  Chem_c <- Chem_c[grep("Homo sapiens", Chem$description), ]

## ---- eval=FALSE--------------------------------------------------------------
#  Chem_c1 <- ms_isoform_resolve(Chem_c)
#  # The next isoform cleaning is optional.
#  # You should and are encouraged to double check and adjust the match table entries.
#  Chem_c2 <- ms_isoform_consolidate(Chem_c1, nread=9, matchtable = "./subfolder/tobe_consolidated.txt")

## ----eval=FALSE---------------------------------------------------------------
#  Chem_c3 <- imprints_rearrange(Chem_c2, nread=9, repthreshold=0.8, countthreshold=1)
#  Chem_c3 <- imprints_rearrange(Chem_c2, nread=9, repthreshold=1, countthreshold=3)

## ----eval=FALSE---------------------------------------------------------------
#  Chem_s <- imprints_normalization(Chem_c3)

## ----eval=FALSE---------------------------------------------------------------
#  Chem_s1 <- imprints_caldiff(Chem_s, reftreatment="G1S")

## ---- eval=FALSE--------------------------------------------------------------
#  Chem_s2 <- imprints_reproducible(Chem_s1)

## ---- eval=FALSE--------------------------------------------------------------
#  Chem_score <- imprints_score(Chem_s1, fdrthreshold=0.05, labelnodes=TRUE, labelcategory="CC+-")
#  Chem_score_mad <- imprints_score(Chem_s1, useMAD=TRUE, nMAD=2.5, labelnodes=TRUE, labelcategory="CC+-")
#  Chem_score_4 <- imprints_score(Chem_s1, format="4", labelnodes=TRUE, labelcategory="CC")
#  
#  # could use the reproducible subset if the previous step was performed
#  Chem_score_r <- imprints_score(Chem_s2, labelnodes=TRUE, labelcategory="CC+-")
#  Chem_score_4r <- imprints_score(Chem_s2, format="4", labelnodes=TRUE, labelcategory="CC")

## ---- eval=FALSE--------------------------------------------------------------
#  # plot out some proteins
#  cyclin <- Chem_s1[grep("[Cc]yclin[- ]",Chem_s1$description),]
#  imprints_barplotting(cyclin, treatmentlevel=c("G1S","S","PM"))
#  # plot out the proteins from one category
#  Chem_CC <- subset(Chem_s1, id %in% subset(Chem_score_4, category=="CC")$id)
#  imprints_barplotting(Chem_CC, treatmentlevel=c("G1S","S","PM"))
#  # plot out all the proteins
#  imprints_barplotting(Chem_s1, treatmentlevel=c("G1S","S","PM"))

## ---- eval=FALSE--------------------------------------------------------------
#  # to first average the data
#  Chem_s1_ave <- imprints_average(Chem_s1)
#  # to check the proteins in PM with Changes to complex
#  Chem_PM_complex <- imprints_complex_mapping(Chem_s1_ave, Chem_score_4,
#                                            treatment="PM", targetcategory=c("NC","CN","CC"),
#                                            organism="Human")
#  # to calculate the correlation of IMPRINTS profiles
#  Chem_PM_complex_corr <- imprints_corr_in_complex(Chem_PM_complex, nread=6)
#  # to remove the redundant complexes with highly similar composition of protein subunits
#  Chem_PM_complex_corr <- imprints_corr_in_complex(Chem_PM_complex, nread=6,
#                                                   removeredundancy=T)

## ---- eval=FALSE--------------------------------------------------------------
#  # It is possible to look at the whole dataset, but it is not recommended to run like below
#  # Chem_s1_corr <- imprints_correlation(Chem_s1_ave, treatmentvector=c("S","PM"))
#  # It is probably more useful to only look at the "Changed" proteins
#  Chem_s1_ave_C <- subset(Chem_s1_ave, id %in% subset(Chem_score_4,
#                                                      category%in%c("CC","NC","CN"))$id)
#  Chem_s1_corr1 <- imprints_correlation(Chem_s1_ave_C, treatmentvector=c("S","PM"))
#  
#  # Then it is possible to plot these "Changed" proteins in a descending order of the correlations, which is controlled by specifying the correlation information
#  Chem_s1_C <- subset(Chem_s1, id %in% subset(Chem_score_4, category%in%c("CC","NC","CN"))$id)
#  imprints_barplotting(Chem_s1_C, treatmentlevel=c("G1S","S","PM"), corrtable=Chem_s1_corr1)

## ---- eval=FALSE--------------------------------------------------------------
#  Chem_PM_de <- imprints_diffExp(Chem_s, contrast="PM-G1S")
#  Chem_S_de <- imprints_diffExp(Chem_s, contrast="S-G1S")

