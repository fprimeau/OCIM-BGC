function [co2,k0,g_k0,g_co2,gg_co2] = eqco2(dic,arg1)
% unpack parameters for co2 system chemistry
alk  = arg1.alk  ;
salt = arg1.salt ;
temp = arg1.temp ;
pres = arg1.pres ;
si   = arg1.si   ;
po4  = arg1.po4  ;

% co2 system
a = CO2SYS(alk,abs(dic),1,2,salt,temp,temp,pres,pres,si,po4,1,4,1);
ax = CO2SYS(alk,abs(dic)+sqrt(-1)*eps^3,1,2,salt,temp,temp,pres,pres,si,po4,1,4,1);
g_k0 = imag(ax(:,8)./ax(:,4))/eps^3;

% concentration of co2 in umol/kg
co2   = a(:,8);
g_co2 = imag(ax(:,8))/eps^3;
% co2 solubility k0 
pco2  = a(:,4); % uatm
k0    = co2./pco2; % mol/kg/atm
%if (nargout>2)
%  % Revelle factor
%  R = a(:,14);
%  g_co2 = R.*co2./dic;
%end
if (nargout>3)
    a = CO2SYS(alk,abs(dic)+sqrt(-1)*eps^3,1,2,salt,temp,temp,pres,pres,si,po4,1,4,1);
    R = a(:,14);
    tmp = R.*co2./dic;
    gg_co2 = imag(tmp)./eps^3;
end
% pH (pH units)
% pH = a(:,3);
