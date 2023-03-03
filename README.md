# IMPRINTS.CETSA
A package primarily designed to process IMPRINTS formats of MS-CETSA result

# What is IMPRINTS.CETSA For?
This package is primarily designed to process IMPRINTS formats of MS-CETSA result, which was first introduced in [Dai et al. Modulation of Protein-Interaction States through the Cell Cycle](https://www.cell.com/cell/fulltext/S0092-8674(18)30397-0).  
[IMPRINTS](https://www.annualreviews.org/doi/10.1146/annurev-biochem-062917-012837) stands for the Integrated Modulation of PRotein INteraction States, is a novel way of sysmematically dissecting the fundamental biochemical modulations during the transitions of cellular states. It is an extension of original [MS-CETSA](https://www.cetsa.org/about). Reseachers who are interested to use MS-CETSA to identify and validate the binding of ligands/small molecules to target proteins on the proteome scale can refer to [mineCETSA package](https://github.com/nkdailingyun/mineCETSA) for data processing.  
We have also make available an accompanying user-friendly shinny-based data visualization package [mineCETSAapp](https://github.com/mgerault/mineCETSAapp).  

The main tasks can be accompolished in this package include:  
1. Read in the quantitative data (Search result txt.file exported from Proteome Discoverer, PD) into R/Rstudio  
2. Clean up data, mainly to remove the entries without quantitative information  
3. Rearrange data into one unique protein per row format  
4. Normalize data based on the principle of minimal variance should be controlled for different samples under each heating condition (i.e.,treatment/replicate combination)  
5. Calculate the relative fold change between treatment for each temperature/replicate combination  
6. Calculate the protein abundance score and stability score, segregate proteins into a few possible categories, such as NN, NC, CN or CC  
7. Visualize the data using customizable barplot layout, for either all the data or the subsets  
8. Map proteins onto known complexes, calculate the correlation of IMPRINTS profiles among the protein complex subunits  
9. Extrapolate the IMPRINTS signal onto a reference melting curve, for further visual and mathematical control  

--------------------------------------------------------------------------------------------

## Prerequisites of installing IMPRINTS.CETSA package
* R version >= 4.0.0, preferably the most updated version
* Rstudio version > 1.0, preferably the most updated version
* Imported packages:  
    arrayQualityMetrics,
    Biobase,
    dplyr (>= 1.0.0),
    drc (>= 3.0),
    fdrtool,
    GGally,
    ggplot2 (>= 3.2.0),
    ggpubr,
    ggrepel (>= 0.8.0),
    graphics,
    grDevices,
    grid,
    gridExtra,
    gridSVG,
    gtools,
    limma,
    MESS,
    methods,
    mice,
    Nozzle.R1,
    plyr (>= 1.8.0),
    RColorBrewer,
    readr (>= 2.0.0),
    reshape2 (>= 1.4.0),
    scales,
    stats,
    tibble (>= 3.0.0),
    tidyr (>= 1.2.0),
    utils,
    VennDiagram,
    VIM,
    XML
During the package installation process from local source, when these dependent packages are not available, complains would rise, use `install.packges()` or `BiocManager::install()` command to manually install the required supporting packages (as specified in the warning message), till there is no more complain.  
