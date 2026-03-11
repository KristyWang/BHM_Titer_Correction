# A modelling framework to improve antibody titer estimation from dilution series data: application to RSV Foci Reduction Neutralization Tests

Yan Wang, Qianli Wang, Chris Wymant, Junyi Zou, Lan Yi, Meng Xu, James A. Hay, Hongjie Yu

## Software and packages
The code was written in R (version 4.4.0). The following packages are needed to run the scripts:<br>
• readxl<br>
• data.table<br>
• dplyr<br>
• tidyverse<br>
• ggplot2<br>
• ggpubr<br>
• patchwork<br>
• RColorBrewer<br>
• lme4<br>
• lmerTest<br>
• MASS<br>
• VGAM<br>
• drc<br>
• rstan<br>
• bayesplot<br>

## A quick start guide
To facilitate adoption of our framework by laboratories and researchers, we provide a quick start guide (https://kristywang.github.io/BHM_Titer_Correction/
) with detailed simulation data. To successfully run the scripts in the guide, you should place the following files in the same folder:<br>
1. “simu_data.xlsx”: Simulated dilution series data mirroring the structure of data from respiratory syncytial virus (RSV) foci reduction neutralization tests (FRNTs).<br>
2. “Stan_Model.R”: Stan model code used to estimate RSV neutralizing antibody titers while correcting for batch effects and an assay-level covariate (i.e., virus working dilution).<br>
3. “Quick_start_guide.Rmd”: An R Markdown document demonstrating how to run the model, process the data, and visualize the resulting estimates.

## Folder structure:
All R scripts used to reproduce the results presented in the manuscript are also provided. The scripts are written assuming the following directory structure:<br>
1.Data<br>
2.Code<br>
3.Output<br>
Where <br>
• All data needed to perform the analyses should be stored in folder '1. Data';<br>
• All the R scripts are reported in folder '2. Code';<br>
• Results produced by the code will be stored in folder '3. Output'.<br>

## Specification for each R script to reproduce the results: 
### “1_Func_Karber_and_4PL.R”:
This R script defines two custom functions for estimating 50% neutralizing antibody (nAb) titers from RSV Foci Reduction Neutralization Test (FRNT) data. The functions include:<br>
• Karber_function – estimates titers using the Kärber formula<br>
• DRM_function – fits a four-parameter logistic (4PL) model to compute titers<br>
Both functions accept a data frame containing dilution series data and calculate antibody titers based on user-specified column names.
### “2_Model_Stan.R”:
This R script defines the Bayesian hierarchical modeling framework described in the manuscript, implemented via two Stan models:<br>
• Model 1 estimates the effect of each working virus stock using only virus control (VC) data. The objective is to isolate the contribution of different virus stocks to the observed foci counts in the VC (i.e., background signal), while accounting for batch-to-batch variability through a random effect.<br>
• Model 2 estimates neutralizing antibody (nAb) titers for individual serum samples. It uses dilution series data from assay controls—including the internal positive control (PC) and the International Standard (IS)—to estimate and correct for batch effects. The model then analyzes dilution series from serum samples to generate accurate nAb titer estimates, appropriately adjusted for batch variability.<br>
### “3_Simu_Func.R”:
This R script defines a custom function to simulate dilution series data. The generated datasets include VC data, PC data, IS data at different concentrations (denoted as IS1 and IS2), and serum sample data.<br>
### “4_Simulation_Study.R”:
This R script contains the complete workflow for conducting the simulation study. It is organized into four main parts:<br>
• Generates dilution series data that replicate key features of real RSV FRNT experimental data.<br>
• Applies standard titer estimation methods, including the Kärber formula and the 4PL model, to the simulated data.<br>
• Fits Bayesian hierarchical modelling framework to the simulated data.<br>
• Generates figures to visualize and interpret the simulation results, facilitating method comparison and performance assessment.<br>
### “5_Clean_experimental_data.R”:
This R script processes and integrates raw experimental data generated from RSV FRNT assays, transforming it into a clean, analysis-ready format compatible with downstream modeling workflows.<br>
### “6_Measurement_uncertainty.R”:
This R script analyzes variability in RSV FRNT measurements using experimental data from assay controls. The objectives are:<br>
• To characterize measurement variability in both the raw experimental readouts (i.e., foci counts) and the derived foci reduction values.<br>
• To quantify the contributions of within-batch and between-batch variability to titer estimates obtained using standard methods, including the Kärber formula and the 4PL model.<br>
### “7_Fit_experimental_data.R”:
This script applies the proposed modeling framework to real RSV FRNT experimental data, alongside standard methods such as the Kärber formula and the 4PL model. It also generates figures to summarize, visualize, and interpret the analysis results.<br>
### “8_Serological_analysis.R”:
This script evaluates how different titer estimation methods impact RSV serological measures at the population level. Analyses include age-stratified seroprevalence, geometric mean titers (GMTs), fold-rises, and seroconversion rates.<br>
