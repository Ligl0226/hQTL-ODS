# hQTL-ODS
Powerful one-dimensional scan to detect heterotic QTL.

---

## Table of Contents

- [Overview](#overview-and-features)
- [System Requirements](#system-requirements)
- [Installation Guide](#installation-guide)
- [Usage and Demo](#usage-and-demo)
- [License](#license)
- [Authors](#authors)
- [Contact](#contact)
- [Citation](#citation)

## Overview and features:
- The heterotic effects were modeled unbiasedly according to the original definition, i.e., the dominance effect of a locus and its epistatic interaction effect with all other loci were considered, independent of their effect sizes. [http://www.nature.com/articles/ng.3974]
- Fast and powerful implementation of the one-dimensional scan to detect heterotic QTL
- Scanning the individual heterotic effect component of the detected hQTL
- Not only for mid-parent heterosis (MPH) but also for better-parent heterosis (BPH)
- Also, support for analyses of the combination of different effects, i.e., additive + cumulative additive-by-additive effect of a locus, et al.

## System Requirements

### Hardware requirements
- No non-standard hardware required
- Runs on a standard computer with enough RAM to support the in-memory operations

### Software requirements
The hQTL-ODS model was implemented as easy-to-use R functions that are run in the R software environment. R can be freely downloaded from http://www.r-project.org. We also recommend the integrated development environment RStudio, which is also freely available at http://www.rstudio.com.
- R (≥ 4.1.1)
- Packages: Rcpp, RcppEigen, gaston

### OS Requirements
This package is supported for Windows and Linux. The package has been tested on the following systems:
- Windows 10/11
- Rocky Linux 9.6 (Blue Onyx)

## Installation Guide
You don't need to install anything.
Load the R scripts into your R environment:
```r
source("./R/CalKinMat.R")
source("./R/MPH_hQTL_ODS.R")
source("./R/hQTL_component_effects.R")
source("./R/BPH_hQTL_ODS.R")
source("./R/QTL_ODS.R")
```

## Usage and Demo
Please refer to [R01.ODS_code_domoo.R](R01.ODS_code_domoo.R) for details. All expected outputs from the demo dataset can be generated within 2 hours on a standard desktop with 10 CPU cores.

## License
hQTL-ODS is distributed under the [GPL-3.0 license](LICENSE.txt).

# Phenotype
The phenotypic data of all parents and hybrids and the simulated data sets collected in the study of "Powerful one-dimensional scan to detect heterotic QTL" are provided in this folder.

## Citation
Powerful one-dimensional scan to detect heterotic QTL

Li G, Schmidt RH, Zhao Y, Reif JC, Jiang Y. Powerful one-dimensional scan to detect heterotic QTL. bioRxiv [Internet]. 2025; Available from: [https://www.biorxiv.org/content/early/2025/05/12/2025.05.06.652159]

## Authors
Guoliang Li and Yong Jiang

## Contact
lig@ipk-gatersleben.de 

