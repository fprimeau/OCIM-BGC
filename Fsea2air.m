function vout = Fsea2air(par, Gtype)
    grd  = par.grd;
    M3d  = par.M3d;
    iwet = par.iwet;
    % get first layer index
    tmp  = M3d;
    tmp(:,:,2:end) = 0;
    isrf = find(tmp(iwet));

    vSST = par.modT(iwet(isrf));
    vSSS = par.modS(iwet(isrf));
    
    tmp = M3d*0;
    tmp(:,:,1) = par.kw;
    kw = tmp(iwet(isrf));
    tmp = M3d*0;
    tmp(:,:,1) = par.P;
    P = tmp(iwet(isrf));

    if strcmp(Gtype,'CO2')
        co2syspar = par.co2syspar;
        
        tmp = M3d+nan;
        tmp(iwet) = par.DIC;
        vDICs = tmp(iwet(isrf));
        
        scco2 = 2073.1 - 125.62*vSST + 3.6276*vSST.^2 - 0.043219*vSST.^3;
        kw    = kw.*sqrt(660./scco2);
        KCO2  = kw/grd.dzt(1);
        % uatm
        pco2atm = par.pco2_air(1); 
        %
        % co2surf unit umol/kg;
        % k0 unit mol/kg/atm
        [co2surf,k0] = eqco2(vDICs,co2syspar);
        co2surf = co2surf;    % umol/kg;
        co2sat  = k0*pco2atm; % umol/kg
        tmp = M3d*0;
        tmp(iwet(isrf)) = KCO2.*(co2sat - co2surf);
        JgDIC = tmp(iwet);
        vout.JgDIC = JgDIC*1024.5/1000; % umole/kg/s to mmol/m^3/s
        
        % Gradient
        [co2surf,k0,g_k0,g_co2] = eqco2(vDICs,co2syspar);
        tmp = M3d*0;
        tmp(iwet(isrf)) = KCO2.*(g_k0.*pco2atm - g_co2)*1024.5/1000;
        vout.KG = d0(tmp(iwet));
        % Hessian
        [co2surf,k0,g_k0,g_co2,gg_co2] = eqco2(vDICs,co2syspar);
        tmp = M3d*0;
        tmp(iwet(isrf)) = -KCO2.*gg_co2*1024.5/1000;
        vout.KGG = tmp(iwet);    
    end
    %
    if strcmp(Gtype,'O2')
        KO2 = M3d*0;
        sco2 = 1638.0 - 81.83*vSST + 1.483*vSST.^2 - 0.008004*vSST.^3;
        kw = kw.*sqrt(660./sco2);
        KO2(iwet(isrf)) = kw/grd.dzt(1);
        vout.KO2 = d0(KO2(iwet));
        
        o2sat = 0*M3d;
        o2sat(iwet(isrf)) = 1000*o2sato(vSST,vSSS).*P; % mmol/m^3
        vout.o2sat = o2sat(iwet);
    end
end

function o2sat = o2sato(T,S)
% Computes the oxygen saturation concentration at 1 atm total
% pressure
% in mol/m^3 given the temperature (t, in deg C) and the salinity
% (s,
% in permil).
%
% FROM GARCIA AND GORDON (1992), LIMNOLOGY and OCEANOGRAPHY.
% THE FORMULA USED IS FROM PAGE 1310, EQUATION (8).
%
% *** NOTE: THE "A3*TS^2" TERM (IN THE PAPER) IS INCORRECT. ***
% *** IT SHOULDN'T BE THERE.                                ***
%
% o2sato IS DEFINED BETWEEN T(freezing) <= T <= 40(deg C) AND
% 0 permil <= S <= 42 permil
%
% CHECK VALUE:  T = 10.0 deg C, S = 35.0 permil,
% o2sato = 0.282015 mol/m^3
%
    A0  = 2.00907;
    A1  = 3.22014;
    A2  = 4.05010;
    A3  = 4.94457;
    A4  = -2.56847e-1;
    A5  =  3.88767;
    B0  = -6.24523e-3;
    B1  = -7.37614e-3;
    B2  = -1.03410e-2;
    B3  = -8.17083e-3;
    C0  = -4.88682e-7;
    TT  = 298.15-T;
    TK  = 273.15+T;
    TS  = log(TT./TK);
    TS2 = TS.^2;
    TS3 = TS.^3;
    TS4 = TS.^4;
    TS5 = TS.^5;
    CO  = A0 + A1*TS + A2*TS2 + A3*TS3 + A4*TS4 + A5*TS5+...
          S.*(B0 + B1*TS + B2*TS2 + B3*TS3)+...
          C0*(S.*S);
    o2sat = exp(CO);
    %
    %  Convert from ml/l to mol/m^3
    %
    o2sat = (o2sat/22391.6)*1000.0;
end

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
    ax = CO2SYS(alk,abs(real(dic))+sqrt(-1)*eps^3,1,2,salt,temp,temp, ...
                pres,pres,si,po4,1,4,1);
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
        a = CO2SYS(alk,abs(real(dic))+sqrt(-1)*eps^3,1,2,salt,temp, ...
                   temp,pres,pres,si,po4,1,4,1);
        R = a(:,14);
        tmp = R.*co2./dic;
        gg_co2 = imag(tmp)./eps^3;
    end
    % pH (pH units)
    % pH = a(:,3);
end