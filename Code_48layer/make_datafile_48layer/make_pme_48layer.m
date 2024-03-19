% make pme from 48layer using auxiliary data from Tim's paper (48layer model)

clc; clear all; close all


addpath('../../DATA/BGC_48layer/');
load OCIM2_CTL_He_48layer.mat           % Convergence & yr-1

XT3d     = output.grid.XT3d;
YT3d     = output.grid.YT3d;
Smod     = output.salt;                 % psu [10-3 g/g] = [mg/g]
Sflux    = output.saltflux;             % positive into the ocean, [mg/m^2/s]
dVt      = output.grid.dVt;                  % [m^3]
dAt      = output.grid.dAt;
msk = output.M3d;
spa = 365*24*60^2;
spd = 24*60^2;
TRdiv = -TR/spa;                        % convergence & yr-1  -----> divergence & s-1 
iwet = find(msk(:));
tau  = spd*30;                          % 30 days to seconds


% surface restoring operator
d0 = @(x) spdiags(x(:),0,length(x(:)),length(x(:)));
L = msk;
L(:,:,2:end) = 0;
L = d0(L(iwet));              

%-----------Sstar, Sbar, and pme ---------------------%
Sstar       = msk + nan;
Sstar(iwet) = tau*(TRdiv+(1/tau)*L)*Smod(iwet);

msk_sfc            = msk;
msk_sfc(:,:,2:end) = 0;
isrf = find(msk_sfc(:));

Sbar = sum(msk_sfc(iwet).*Smod(iwet).*dVt(iwet))/sum(msk_sfc(iwet).*dVt(iwet));   
pme  = ((Smod-Sstar)/Sbar).*msk_sfc.*(1/tau);            
pme_integral  = sum(dAt(isrf) .* pme(isrf));

% Surface integral -----> not exact to be zero. 
% Make integral of pme close to zero.

pme_new          = pme;
constant_pme     = sum(dAt(isrf).*pme(isrf))/sum(msk_sfc(isrf).*dAt(isrf));
pme_new(isrf)    = pme(isrf) - constant_pme;
pme_new_integral = sum(dAt(isrf).*pme_new(isrf));
                                            

%Save pme file in DATA directory
fileName = 'pme_91x180_48layer.mat'
directory = '../../DATA/BGC_48layer'
filePath = fullfile(directory, fileName);
save(filePath, 'pme_new');
