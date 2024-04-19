# Carbon-Export
The biogeochemical model that couples the cycling of phosphorus, carbon, and oxygen is used to obtain global scale organic carbon export flux. 

driver.m is the main driver

neglogpost.m computes the negative of the logarithm of the posterior pdf for the parameters conditioned on the DIP, DIC, ALK, DOC, and O2 data.

uptake_C.m computes the rate of uptake of DIC given the DIP uptake as well as the 1st and 2nd derivatives w.r.t. the model parameters

eqPcycle.m computes the equilibrium state of the P-cycle model as well as the 1st and 2nd derivaties of the equilibrium solution w.r.t. the model parameters

eqCcycle.m computes the equilibrium of the C-cycle model state as well as the 1st and 2nd derivative of the solution w.r.t the model parameters

eqOcycle.m computes the equilibrium of the O-cycle model state as well as the 1st and 2nd derivative of the solution w.r.t the model parameters

Fsea2air.m computes the air-sea gas exchange for CO2 and O2.

PackPar.m Organize the tunable parameters and assign each parameter an index

PrintPar.m Print the parameters to a log file at each iteration

ResetPar.m Reset a parameter if the optmization routine suggests it a strange value. This function is only used in the first 10 iteration, after that the choice of parameters is solely based on the data.

SetPar.m Assign initial values to parameters

SetUp.m Setup the environment and load supporting data.

buildPFD.m build the particle flux divergence operator.

build_co2syspar.m build the parameters for the co2 system.

mkO2P build the oxygen to phosphorus uptake ratio.

mkPIC2P.m build the PIC to POC production ratio.

unitlity scripts:

nsnew.m C.T. Kelly's Newton Armijo solver

mfactor.m Timothy A. Davis' LINFACTOR VERSION 1.1.0, Nov 1, 2007 Copyright 2007, Timothy A. Davis, University of Florida

d0.m makes a sparse diagonal matrix given a 3d field

CO2SYS.m Lewis, E., and D. W. R. Wallace. 1998. Program Developed for CO2 System Calculations. ORNL/CDIAC-105. Carbon Dioxide Information Analysis Center, Oak Ridge National Laboratory, U.S. Department of Energy, Oak Ridge, Tennessee.


