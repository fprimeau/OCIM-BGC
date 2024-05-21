%set_figoptions
% set default figure options. 
% define colors for plotting

%% ------ set up figures ---------

set(groot,'defaultAxesFontName','Times',...
    'defaultAxesFontSize',9,...
    'defaultAxesTickLabelInterpreter','latex',...
    'defaultAxesXMinorTick','on',...
    'defaultAxesYMinorTick','on');
% TEXT PROPERTIES
set(groot,'defaultTextFontName','Times',...
    'defaultTextInterpreter','latex');

colors.maroon 		= [128/255 0 0];
colors.tomato 		= [255/255 99/255 71/255];	% light red-orange
colors.indianred 	= [205/255 92/255 92/255]; % light red-brown
colors.limegreen 	= [50/255 205/255 50/255];
colors.darkgreen 	= [0 100/255 0];
colors.teal 		= [0 128/255 128/255];
colors.aqua 		= [0.2 0.8 0.8];
colors.lblue 		= [0 191/255 255/255];
colors.navy 		= [ 0 0 128/255];
colors.darkmagenta 	= [139/255 0 139/255];
colors.darkgrey 	= [0.3 0.3 0.3];
colors.lgrey 	    = [0.8 0.8 0.8];


%% define 1 line color for each model for comparison line plots
modelcolors.Tz      = [180/255 0 180/255]; % = colors.darkmagenta * 1.295;
modelcolors.GM15    = colors.lblue;
modelcolors.Cell    = colors.darkgreen;
modelcolors.Const   = colors.darkgrey;
