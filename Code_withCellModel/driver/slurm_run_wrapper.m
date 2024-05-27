% slurm_run_wrapper.m
% running this script in the slurm matlab call instead of running driver.m
% directly should suppress the command line prompts from appearing in the
% output file
diary logs/optPCO_constC2P_noGH.out
fprintf('Test optimizing without computing gradient and hessian \n\n')
run driver_test_noGH.m
diary off