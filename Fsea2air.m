function [JgDIC,KG,KGG] = Fsea2air(par, parm, DIC)
grd  = parm.grd;
M3d  = parm.M3d;
iwet = parm.iwet;
% get first layer index
tmp  = M3d;
tmp(:,:,2:end) = 0;
isrf = find(tmp(iwet));
co2syspar = parm.co2syspar;

tmp = M3d+nan;
tmp(iwet) = DIC;
vDICs = tmp(iwet(isrf));
vSST  = parm.modT(iwet(isrf));

[kw,P] = kwco2(M3d,grd);
tmp = M3d*0;
tmp(:,:,1) = kw;
kw = tmp(iwet(isrf));

scco2 = 2073.1 - 125.62*vSST + 3.6276*vSST.^2 - 0.043219*vSST.^3;
kw    = kw.*sqrt(660./scco2);
KCO2  = kw/grd.dzt(1);
% uatm
pco2atm = parm.pco2_air(1); 
%
% co2surf unit umol/kg;
% k0 unit mol/kg/atm
[co2surf,k0] = eqco2(vDICs,co2syspar);
co2surf = co2surf;    % umol/kg;
co2sat  = k0*pco2atm; % umol/kg
tmp = M3d*0;
tmp(iwet(isrf)) = KCO2.*(co2sat - co2surf);
JgDIC = tmp(iwet);
JgDIC = JgDIC*1024.5/1000; % umole/kg/s to mmol/m^3/s

%
if (nargout > 1)
    [co2surf,k0,g_k0,g_co2] = eqco2(vDICs,co2syspar);
    %
    tmp = M3d*0;
    tmp(iwet(isrf)) = KCO2.*(g_k0.*pco2atm - g_co2)*1024.5/1000;
    KG = d0(tmp(iwet));
end

if (nargout > 2)
    [co2surf,k0,g_k0,g_co2,gg_co2] = eqco2(vDICs,co2syspar);
    %
    tmp = M3d*0;
    tmp(iwet(isrf)) = -KCO2.*gg_co2*1024.5/1000;
    KGG = tmp(iwet);    
end
