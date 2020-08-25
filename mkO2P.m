function vout = mkO2P(par)
% parameters
slopeo   = par.slopeo;
interpo  = par.interpo;
%
iwet = par.iwet;
nwet = par.nwet;
sal  = par.modS;
modT = par.modT;
smsk = par.M3d;
smsk(:,:,2:end) = 0;
isrf = find(smsk(iwet));
dVs = par.dVt(iwet(isrf));
surface_mean = @(x) sum(x(isrf).*dVs)/sum(dVs);

% compute the mean of the regressor variable
Z  = modT(iwet);
mu = surface_mean(Z);
Delta = sqrt(surface_mean((Z-mu).^2));

% standardize the regressor variables
ZR = (Z-mu)/Delta; par.ZR = ZR;
%
vout.O2P = slopeo*ZR + interpo;
vout.dO2Pdslopeo = ZR;
vout.dO2Pdinterpo = sparse(nwet,1)+1;