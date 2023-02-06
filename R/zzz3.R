.onLoad <- function(...) {

  message <- paste("Welcome to IMPRINTS.CETSA! This is a sister package of mineCETSA!",
    "Now it's the time to start mining your IMPRINTS-CETSA data!",
    "The vignette below explain in detail how to use this package",
    "Access the respective vignettes by key in: ",
    "vignette('IMPRINTS-CETSA-vignette', package='IMPRINTS.CETSA')", sep="\n")

  packageStartupMessage(message)
}
