% slurm_run_wrapper.m
% running this script in the slurm matlab call instead of running driver.m
% directly should suppress the command line prompts from appearing in the
% output file
diary logs/optPCO_GM15_N23in_testNPP.out
%fprintf('----------------------------------- \n\n')
dbstop if error
date_start = datetime('now');
fprintf('Start Date: %s\n',date_start)
fprintf('Reoptimizing phosphate-dependent C:P model using the SetUp input fields from the Nature paper (Wang et al., 2023) using initial guess on CO and params from from hojongs previous run. smoothed DIP_obs and Temp obs; No gm15 satellite npp conversions (instead constant p2c = 1/117). Testing if satellite NPP conversion to Pnpp causes issues with optimization. \n\n')
run driver_SetUpv2.m
date_end = datetime('now'); fprintf('Complete Date: %s\n',date_end)
fprintf('Total elapsed time was %s\n',date_end - date_start)
diary off