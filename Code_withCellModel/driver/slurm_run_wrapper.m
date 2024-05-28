% slurm_run_wrapper.m
% running this script in the slurm matlab call instead of running driver.m
% directly should suppress the command line prompts from appearing in the
% output file
diary logs/GHtest_sigC_rO2C.out
fprintf('Test fixed rO2C gradient \n\n')
run driver_GHtest_tmp.m
diary off