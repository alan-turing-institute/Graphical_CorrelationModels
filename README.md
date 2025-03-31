#  A graphical framework for interpretable correlation matrix models



This repository contains data and code to compute models for correlation matrices with a user-defined graphical structure. The graphical structure makes correlation matrices interpretable and avoids the quadratic increase of parameters as a function of the dimension. We suggest an automatic approach to define a prior using a natural sequence of simpler models within the Penalized Complexity framework for the unknown parameters in these models. 

## Description

This folder contains data and code to accompany the paper on 'A Graphical Framework for Interpretable
Correlation Matrix Models for Multivariate
Regression'. 

We have included R  scripts to compute the three applications presented in the paper. 

	

##  Repository Structure
```
Code/
├── GermanyData.Rdata  # data for Multivariate Disease Mapping
├── Simple_Model_Biomarkers.R  #  Biomarkers data generated - original data not disclosed.
├── Multivariate_ Disease_Example.R # Multivariate disease mapping 
├── RGEN.R #  pc priors functions 
└── mbesag.R # auxiliary functions

Code/Simulations/ 

├── corGraphs.tar.gz  #  R package with pc priors function code in C to improve computational times 
├── Define_lambda.R  # code for Fig. 3 and Fig.4 
├── 2Y/ # code for all the simulations for the different sample sizes for section 5.1.2 
├── 4Y/ #  code for all the simulations for the different sample sizes for section 5.1.3
├── LowCor #  code for all the simulations with low correlations reported in the appendix
├── results_2Y.R  # code for table 2
└── results_4Y.R  #  code for tables 3 and 4


```

## Getting Started

### Dependencies

* The code runs on R 4.4.1 (2024-06-14)
* INLA_24.10.05 

### Executing program

* Run the files by file order as reported. 


## R scripts Description 

*GermanyData.Rdata: is the dataset that contains expected and observed data for 4 cancer cases observed in Germany	

*Simple_Model_Biomarkers.R: it generates four variables representing the epigenetic clocks and fits a simple multivariate normal model with no covariates.

*Multivariate_ Disease_Example.R: it computes the five bym models and the codes for the maps.

*RGEN.R: contains specific functions to fit the correlation models and the code to generate the graphs, to source only.

*mbesag.R: contains function specific to the multivariate disease application, to source only.

*corGraphs.tar.gz: this is the draft of the R package where the code from RGEN.R has been implemented in C; it is still under development. 

*Define_lambda.R:  this code produces the simulations for different levels of lambda and produces the plots reported in     Fig. 3 and Fig.4 

*2Y/: this folder contains all the code for all the simulations for the different sample sizes for section 5.1.2; it contains the simulations for 30 and 1000 subjects, different lambda levels and misspecified and unspecified correlation structures. 

*4Y/: this folder contains all the code for all the simulations for the different sample sizes for section 5.1.3; it contains the simulations for 30 and 1000 subjects, different lambda levels and unspecified correlation structures.

*results_2Y.R: this code creates the table to summarize the results in table 2.

*results_4Y.R: this code creates tables to summarize the results in tables 3 and 4.  

*LowCor/: contains the simulations with lower correlation values  


## Authors


[Anna Freni Sterrantino](mailto:afrenisterrantino@turing.ac.uk).


## License

This work is licensed under the MIT license (code) and Creative Commons Attribution 4.0 International license (for documentation). You are free to share and adapt the material for any purpose, even commercially, as long as you provide attribution (give appropriate credit, provide a link to the license, and indicate if changes were made) in any reasonable manner, but not in any way that suggests the licensor endorses you or your use, and with no additional restrictions.


