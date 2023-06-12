# IMPRINTS.CETSA
A package primarily designed to process the IMPRINTS format of MS-CETSA results

# What is IMPRINTS.CETSA For?
This package is primarily designed to process the IMPRINTS format of MS-CETSA results, first introduced in [Dai et al. 2018 Cell](https://www.cell.com/cell/fulltext/S0092-8674(18)30397-0).  

[IMPRINTS](https://www.annualreviews.org/doi/10.1146/annurev-biochem-062917-012837), which stands for the Integrated Modulation of PRotein INteraction States, is a novel way to systematically analyze the underlying biochemical modulations during cellular state transitions. It is an extension of the original [MS-CETSA](https://www.cetsa.org/about). Researchers interested in using MS-CETSA to identify and validate the binding of ligands/small molecules to target proteins on a proteome scale can refer to the [mineCETSA package](https://github.com/nkdailingyun/mineCETSA) for data processing.  

We have also provided a companion user-friendly Shiny-based data visualization package [IMPRINTS.CETSA.app](https://github.com/mgerault/IMPRINTS.CETSA.app).  

The main tasks that can be accomplished in this package include:  
1. Import the quantitative data (search result txt.file exported from Proteome Discoverer, PD) into R/Rstudio  
2. Clean up the data, mainly to remove the entries without quantitative information  
3. Rearrange data into one unique protein per row format  
4. Normalize data based on the principle of minimal variance should be controlled for different samples under each heating condition (i.e.,treatment/replicate combination)  
5. Calculate the relative fold change between treatments for each temperature/replicate combination  
6. Calculate the protein abundance score and stability score, segregate proteins into a few possible categories, such as NN, NC, CN or CC  
7. Visualize the data using a customizable bar plot layout, either for all data or for the subsets  
8. Map proteins to known complexes, calculate the correlation of IMPRINTS profiles among the protein complex subunits  
9. Extrapolate the IMPRINTS signal to a reference melting curve, for further visual and mathematical control  

--------------------------------------------------------------------------------------------

## How to install IMPRINTS.CETSA ?  
First go to Rstudio and check that you have R >= 4.0
Moreover, you'll need to install two packages from Bioconductor with this commands:

```c
if(!requireNamespace("BiocManager", quietly = TRUE)){
   install.packages("BiocManager") 
}
BiocManager::install(c("limma", "arrayQualityMetrics"))
```

When all of this is done, type and run the following commands in R console:

```c
if(!requireNamespace("devtools", quietly = TRUE)){
   install.packages("devtools")
} 
devtools::install_github("nkdailingyun/IMPRINTS.CETSA")
```

IMPRINTS.CETSA is now installed and you can load it with:

```
library(IMPRINTS.CETSA)
```
