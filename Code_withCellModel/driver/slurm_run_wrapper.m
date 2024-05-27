% slurm_run_wrapper.m
% running this script in the slurm matlab call instead of running driver.m
% directly should suppress the command line prompts from appearing in the
% output file
diary logs/optPCOCell_v1continued.out
fprintf('Continuing optPCOCell_v1 optimization after TIMEOUT \n\n')
run driver.m
diary off