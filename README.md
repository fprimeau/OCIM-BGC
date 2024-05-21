# OCIM-BGC-48layer Model

## Project Status
*NOTE*: This code is a *work in progress* and is not intended for use outside of the research group. 

## Contents
### Code_48layer: src, driver, and make_datafile
- `src`: Source code for BGC modeling
- `driver`: Main driver for BGC modeling
- `make_datafile`: Make the datafile for 48layer model (i.e., pme, NPP from satellite data...)

### Code_withCellModel: src, driver
- `src`: Source code for BGC modeling; includes mechanistic stoichiometry model; 24 vertical layers
- `driver`: Main driver scripts to run the BGC model

### DATA: Supporting data used to run the inverse model
Some large data files (e.g., climatological mean NPP data and Transport operator of 48layer) are not uploaded to GitHub. Please refer to the DFS-L in greenplanet. We will also move all the raw data and supporting data files to another drive.

## References
The file `atmhist.txt` is from the supplementary materials of the following article:
```
@article{JGRC:JGRC12521,
author = {Graven, H. D. and Gruber, N. and Key, R. and Khatiwala, S. and Giraud, X.},
title = {Changing controls on oceanic radiocarbon: New insights on shallow-to-deep ocean exchange and anthropogenic CO2 uptake},
journal = {Journal of Geophysical Research: Oceans},
volume = {117},
number = {C10},
issn = {2156-2202},
url = {http://dx.doi.org/10.1029/2012JC008074},
doi = {10.1029/2012JC008074},
pages = {n/a--n/a},
keywords = {Oceans, Numerical modeling, Carbon cycling, Chemical tracers, Gases, anthropogenic CO2, bomb radiocarbon, ocean modeling, transient tracers},
year = {2012},
note = {C10005},
}
```

## To Do List
1. Sensitivity test for the vertical structure of NPP.
2. Tune the piston velocity and reoptimize the BGC model parameters to find DIC in the preindustrial period.
3. Clean up the code and merge it with the 24layer model.
