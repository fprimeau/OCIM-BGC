% slurm_run_wrapper.m
% running this script in the slurm matlab call instead of running driver.m
% directly should suppress the command line prompts from appearing in the
% output file
diary logs/optPCO_Tz_v1.out
fprintf('C2P as a linear function of normalized temperature \n\n')
run driver.m
diary off