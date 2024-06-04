% slurm_run_wrapper.m
% running this script in the slurm matlab call instead of running driver.m
% directly should suppress the command line prompts from appearing in the
% output file
diary logs/optPCO_Tz_v2.out
fprintf('----------------------------------- \n\n')
datetime('now')
fprintf('continue optim after fixing gradient for ccT and ddT \n\n')
run driver.m
diary off