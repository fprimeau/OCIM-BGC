% slurm_run_wrapper.m
% running this script in the slurm matlab call instead of running driver.m
% directly should suppress the command line prompts from appearing in the
% output file
diary logs/optPCO_GM15_v1continued.out
fprintf('Continue optPCO_GM15_v1 after slurm Timeout \n\n')
run driver.m
diary off