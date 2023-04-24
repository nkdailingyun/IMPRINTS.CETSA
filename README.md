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
