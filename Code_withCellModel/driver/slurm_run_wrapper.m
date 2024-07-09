% slurm_run_wrapper.m
% running this script in the slurm matlab call instead of running driver.m
% directly should suppress the command line prompts from appearing in the
% output file
clear all; close all; clc;
diary logs/optPCO_GM15_N23in_testNPPp2c_Nature_ccdd.out
%fprintf('----------------------------------- \n\n')
dbstop if error
date_start = datetime('now');
fprintf('Start Date: %s\n',date_start)
fprintf('Reoptimizing phosphate-dependent C:P model using the SetUp input fields from the Nature paper (Wang et al., 2023) using initial guess on CO and params from from hojongs previous run. no smoothing for Temp obs or DIP_obs; use optimized p2c cc and dd parameters from Nature 2023 for the satellite npp conversions. Testing if cc and dd change for slightly differnt satellite p2c conversion. using p2c consistent with expected result. \n\n')
run driver_SetUpv2.m
date_end = datetime('now'); fprintf('Complete Date: %s\n',date_end)
fprintf('Total elapsed time was %s\n',date_end - date_start)
diary off