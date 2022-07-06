%% test_alphaPhoto
% compare photosynthesis rate from CellPhoto vs in-line funciton alphaPhoto 
% in CellCNP (chebychev polynomial approximation)

%% add path to OCIM-BGCv2
addpath('/Users/megansullivan/Documents/UC Irvine/GitHub/OCIM-BGC-Cell')

%% load light input field
load('Cellmodel_input_obs.mat')

%spd  = 24*60^2 ;  	% seconds per day
%spa  = 365*spd ;  	% seconds per year
%rho   = 1024.5       ; % seawater density [kg/m^3];
%permil    = rho*1e-3 ; % density*[1mmol/10^3 umol]; used for conversion from umol/kg to mmol/m3;

grd.xt = obs.gridx ;
grd.yt = obs.gridy ;

% load in only the top layers
M3d_EZ = obs.wetpoints(:,:,1:2);

%par.po4obs = obs.po4(:,:,1:2);
%par.no3obs = obs.no3(:,:,1:2) ;
%par.Temp = obs.temp(:,:,1:2) ;
par.PARobs = obs.PAR(:,:,1:2) ;

iprod = find(M3d_EZ); %production in top two layers
%P0 = par.po4obs(iprod)./10^6;   %  convert [ mmol/m^3 --> mol/L]
%N0 = par.no3obs(iprod)./10^6;   %  convert [ mmol/m^3 --> mol/L]
%T0 = par.Temp(iprod);            % [units: degC]
Irr0 = par.PARobs(iprod);        % [units: umol photon m^-2 s^-1]

%% first test mean value
Irravg = mean(par.PARobs(iprod));
fprintf('Mean PAR is %6.2f umol photons m^-2 s^-1 \n',mean(par.PARobs(iprod)))

% chebychev polynomial approximation
% alphaPhoto function in CellCNP

% rational Chebyshev
        I0 = Irravg/200.0; % scale by roughly magnitude of irradiance (use avg irradience of box, assuming exponential decay)
        LI = 2.0;   % Length of
        ITrans = (I0-LI)./(I0+LI); % transforms interval to -1 to 1 . interpolation works on this variabble

        Ch0 = 0.19590903;
        Ch1 = 0.10462839;
        Ch2 = -0.05337273;
        Ch3 = 0.01467971;

        aPhotoavg_approx = Ch0+Ch1*(ITrans)+Ch2*(2*ITrans.^2-1)+Ch3*(4*ITrans.^3-3*ITrans);

fprintf('Photosynthetic rate at 25C is %6.4f (using chebyshev polynomial approximation) \n', aPhotoavg_approx)

% run function CellPhoto
[F1opt,fphoto_star] = CellPhoto(Irravg);
aPhotoavg = fphoto_star; %photosynthetic rate, temperature independent
fprintf('Photosynthetic rate at 25C is %6.4f (using CellPhoto function) \n \n', aPhotoavg)


%% -------- for All surface values -----------------------
Irr0 = par.PARobs(iprod);        % [units: umol photon m^-2 s^-1]

% rational Chebyshev polynomial approximation
        I0 = Irr0./200.0; % scale by roughly magnitude of irradiance (use avg irradience of box, assuming exponential decay)
        LI = 2.0;   % Length of
        ITrans = (I0-LI)./(I0+LI); % transforms interval to -1 to 1 . interpolation works on this variabble

        Ch0 = 0.19590903;
        Ch1 = 0.10462839;
        Ch2 = -0.05337273;
        Ch3 = 0.01467971;

        aPhoto_approx = Ch0+Ch1*(ITrans)+Ch2*(2*ITrans.^2-1)+Ch3*(4*ITrans.^3-3*ITrans);

% run function CellPhoto
aPhoto = NaN(size(Irr0));
for ii = 1: length(Irr0)
    [F1opt,fphoto_star] = CellPhoto(Irr0(ii));
    aPhoto(ii) = fphoto_star; %photosynthetic rate, temperature independent
end

%% plot aPhoto_approx vs aPhoto
figure;
plot(aPhoto,aPhoto_approx,'*b'); hold on
plot([0 max(aPhoto)],[0 max(aPhoto)],'-k','Linewidth',2)
xlabel('CellPhoto Function')
ylabel('Chebyshev Polynomial approximation')
title('alphaPhoto (temperature independent photosynthetic rate per unit of the photosynthetic apparatus)')
grid on

%2 subplots. 1) map of aI_approx 2) map of aI

%% plot aPhoto as function of Irr0
figure;
plot(Irr0,aPhoto_approx,'.r'); hold on
plot(Irr0,aPhoto,'.b')
title('CellPhoto photosynthesis rate vs Irradiance')
xlabel('PAR [umol photon m^-^2 s^-^1]')
ylabel('photosynthetic rate per unit photosynthetic apparatus')
grid on
legend('Chebyshev approx','CellPhoto function','location','best')

%% where is aPhoto < 0.05 geographically?

aPhoto_sm = NaN(size(aPhoto));
aPhoto_sm(aPhoto<0.05) = aPhoto(aPhoto<0.05);

aPhoto_smgrid = NaN(size(M3d_EZ));
aPhoto_smgrid(iprod) = aPhoto_sm;

figure;
subplot(2,1,1)
contourf(grd.xt,grd.yt,aPhoto_smgrid(:,:,1)); c1 = colorbar;
title('Surface')

subplot(2,1,2)
contourf(grd.xt,grd.yt,aPhoto_smgrid(:,:,2)); c2 = colorbar;
title('Lower EZ')

% at the surface, alphaPhoto is only less than 0.05 in the Weddell Sea
% in the lower Euphotic layer, alphaPhoto is also below 0.05 in the arctic

