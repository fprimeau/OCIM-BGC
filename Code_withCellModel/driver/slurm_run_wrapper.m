% slurm_run_wrapper.m
% running this script in the slurm matlab call instead of running driver.m
% directly should suppress the command line prompts from appearing in the
% output file
diary logs/optPCO_constC2P_v3.out
fprintf('restart optPCO with new behavior of ResetPar in neglogpost \n\n')
run driver.m
diary off