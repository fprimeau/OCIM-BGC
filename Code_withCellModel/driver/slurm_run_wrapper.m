% slurm_run_wrapper.m
% running this script in the slurm matlab call instead of running driver.m
% directly should suppress the command line prompts from appearing in the
% output file
diary logs/reoptNature_GM15_PCO.out
%fprintf('----------------------------------- \n\n')
date_start = datetime('now');
fprintf('Start Date: %s\n',date_start)
fprintf('Reoptimizing the Nature paper (Wang et al., 2023) using parameters from hojongs previous run and only BGC_2023Nature input data \n\n')
run driver_reoptNature.m
date_end = datetime('now'); fprintf('Complete Date: %s\n',date_end)
fprintf('Total elapsed time was %s\n',date_end - date_start)
diary off