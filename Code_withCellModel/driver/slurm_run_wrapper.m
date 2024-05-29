% slurm_run_wrapper.m
% running this script in the slurm matlab call instead of running driver.m
% directly should suppress the command line prompts from appearing in the
% output file
diary logs/optPCO_Tz_v2.out
fprintf('Restart optPCO_Tz with new ResetPar behavior \n\n')
run driver.m
diary off