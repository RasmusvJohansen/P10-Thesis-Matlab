# P10-Thesis-Matlab

# Table of contents
1. [Requirments](#requirements)
1. [Repository structure](#repository-structure)

## Requirements
This section list of the requirements to run the repository: 
* MATLAB min 2023b
* [YALMIP](https://yalmip.github.io/)
* [MOSEK solver](https://www.mosek.com/)



## Repository structure
This section describes the structure of the repository. The repository contains all the files for running the examinations and synthesis of the project.  All the individual examination and synthesis files corresponds to a specific test in the rapport and can therefore be used to go though each section of the rapport individually. When running the examinations and synthesis images shown in the rapport are generated again and saved in the folder Images. 

The Repository contains the folders `+control, +model and +util`. These folders are used to the respective functions for calculating the control related aspects of the project, and model contains functions for calculating model related aspects of the project. Lastly the util folder contains a mixture of utility functions for converting between units and simulate steps. The functions of the folders can be accessed in a matlab script by calling the folder.functionname(parameters) e.g.
`util.Celsius2Kelvin(20)` would run the function converting celsius to kelvin, stored in the util folder. 

Lastly it should be noted that the program is based around the struct "param" defined in the model.parameters this is used to store calculated parameters of the system, E.g model parameters, weight functions, controller values. And can therefore also be used in functions which requires the system parameters to run. A example of the usage of param can be see in `Examination_SOF_Performance.mlx`