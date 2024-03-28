#  A graphical framework for interpretable correlation matrix models


This repository contains data and code to compute models for correlation matrices with a user-defined graphical structure. The graphical structure makes correlation matrices interpretable and avoids the quadratic increase of parameters as a function of the dimension. We suggest an automatic approach to define a prior using a natural sequence of simpler models within the Penalized Complexity framework for the unknown parameters in these models. 

## Description

This folder contains data and code to accompany the paper on ['A graphical framework for interpretable correlation matrix models']( https://arxiv.org/abs/2312.06289). 

We have included R  scripts to compute the three applications presented in the paper. 



		

##  Repository Structure
```
Code/
├── GermanyData.Rdata  # data for Multivariate Disease Mapping
|
├── Simple_Model_Biomarkers.R  #  biomarkers data generated - original data not disclosed.
├── Multivariate_ Disease_Example.R
├── Longitudinal_Example_Linear.R
├── Longitudinal_Example_nonLinear.R
├── Lambda_Simulations.R
├── RGEN.R  
└── mbesag.R
```

## Getting Started

### Dependencies

* The code runs on R 4.3.2 (2023-10-31)
* INLA_23.09.09 

### Executing program

* Run the files by application. 

## R scripts Description 

*GermanyData.Rdata: is the dataset that contains expected and observed data for three cancer cases observed in Germany	
*Simple_Model_Biomarkers.R: it generates four variables to represent the epigenetic clocks and fits a simple multivariate normal model with no covariates.
*Multivariate_ Disease_Example.R: it computes the five bym models and the codes for the maps.
*Longitudinal_Example_Linear.R: Provides the simulated data and graphs definition for the longitudinal example with two biomarkers and linear time effects.
*Longitudinal_Example_nonLinear.R: provides the simulated data for the longitudinal data with quadratic time effects for two correlated biomarkers.
*Lambda_Simulations: provide the simulation codes for the lambda definition and the code for generating the figures in the paper.
*RGEN.R: contains specific functions to fit the correlation models and the code to generate the graphs to source only.
*mbesag.R: contains function specific to the multivariate disease application, to  source only.


## Authors

[Anna Freni Sterrantino](mailto:afrenisterrantino@turing.ac.uk).


## License

This work is licensed under the MIT license (code) and Creative Commons Attribution 4.0 International license (for documentation). You are free to share and adapt the material for any purpose, even commercially, as long as you provide attribution (give appropriate credit, provide a link to the license, and indicate if changes were made) in any reasonable manner, but not in any way that suggests the licensor endorses you or your use, and with no additional restrictions.


