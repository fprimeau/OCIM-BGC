% slurm_run_wrapper.m
% running this script in the slurm matlab call instead of running driver.m
% directly should suppress the command line prompts from appearing in the
% output file
diary logs/GHtest_sigC_rO2C_O2CT.out
fprintf('testing fixed Hessian for rO2C \n\n')
run driver_GHtest_tmp.m
diary off