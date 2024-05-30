%load_modeloutput_opt
%% load in optimized model output from all runs

    output_dir = '../output/';

    % load Constant model output
	fname = strcat(output_dir,'optPCO_constC2P_v3_CTL_He_PCO_DOC1_DOP0.mat');
	load(fname);
	model_Const = data;
	fxhat = strcat(output_dir,'optPCO_constC2P_v3_CTL_He_PCO_DOC1_DOP0_xhat.mat');
	load(fxhat);
	xhat_Const = xhat;
	clear data xhat

	% load Temperature model output
	fname = strcat(output_dir,'optPCO_Tz_v2_CTL_He_PCO_DOC1_DOP0.mat');
	load(fname);
	model_Tz = data;
	fxhat = strcat(output_dir,'optPCO_Tz_v2_CTL_He_PCO_DOC1_DOP0_xhat.mat');
	load(fxhat);
	xhat_Tz = xhat;
	clear data xhat

	% load GM15 model output
	fname = strcat(output_dir,'optPCO_GM15_CTL_He_PCOv1_DOC1_DOP0.mat');
	load(fname);
	model_GM15 = data;
	fxhat = strcat(output_dir,'optPCO_GM15_CTL_He_PCOv1_DOC1_DOP0_xhat.mat');
	load(fxhat);
	xhat_GM15 = xhat;
	clear data xhat

	% load cell model output fields
	fname = strcat(output_dir,'optPCOCell_CTL_He_PCOCellv1_DOC1_DOP0.mat');
	load(fname);
	model_Cell = data;
	% load optimal parameter values
	fxhat = strcat(output_dir,'optPCOCell_CTL_He_PCOCellv1_DOC1_DOP0_xhat.mat');
	load(fxhat);
	xhat_Cell = xhat;
	clear data xhat
