function vout = mkO2P(parm)
% parameters
slopeo   = parm.slopeo;
interpo  = parm.interpo;
%
iwet = parm.iwet;
nwet = parm.nwet;
sal  = parm.ss;
sst  = parm.sst;
smsk = parm.M3d;
smsk(:,:,2:end) = 0;
isrf = find(smsk(iwet));
dVs = parm.dVt(iwet(isrf));
surface_mean = @(x) sum(x(isrf).*dVs)/sum(dVs);

% compute the mean of the regressor variable
Z = sst(iwet);
mu = surface_mean(Z);
Delta = sqrt(surface_mean((Z-mu).^2));

% standardize the regressor variables
ZR = (Z-mu)/Delta; parm.ZR = ZR;
%
vout.O2P = slopeo*ZR + interpo;
vout.dO2Pdslopeo = ZR;
vout.dO2Pdinterpo = sparse(nwet,1)+1;