function vout = Fsea2air(par, Gtype)
    grd  = par.grd  ;
    M3d  = par.M3d  ;
    iwet = par.iwet ;
    % get first layer index
    tmp  = M3d      ;
    tmp(:,:,2:end) = 0     ;
    isrf = find(tmp(iwet)) ;

    vSST = par.Temp(iwet(isrf)) ;
    vSSS = par.Salt(iwet(isrf)) ;
    
    tmp = M3d*0           ;
    tmp(:,:,1) = par.kw   ;
    kw  = tmp(iwet(isrf)) ;
    tmp = M3d*0           ;
    tmp(:,:,1) = par.P    ;
    P   = tmp(iwet(isrf)) ;

    if strcmp(Gtype, 'CO2')
        
        pco2atm   = par.pco2atm      ;  % uatm
        vDICs     = par.DIC(isrf)    ;
        vALKs     = par.ALK(isrf)    ;
        co2syspar = par.co2syspar    ;
        scco2 = 2073.1 - 125.62*vSST + 3.6276*vSST.^2 - 0.043219*vSST.^3;
        kw    = kw.*sqrt(660./scco2) ;
        KCO2  = kw/grd.dzt(1)        ;
        %
        % co2surf unit umol/kg;
        % k0 unit mol/kg/atm
        [co2surf,k0] = eqco2(vDICs,vALKs,co2syspar) ;
        co2sat = k0*pco2atm    ; %mol/kg/atm -> umol/kg
        tmp    = M3d*0         ;
        tmp(iwet(isrf)) = KCO2.*(co2sat - co2surf)*par.permil ;
        vout.JgDIC = tmp(iwet) ; % umole/kg/s to mmol/m^3/s
        
        % Gradient
        [co2surf,k0,Gout] = eqco2(vDICs,vALKs,co2syspar) ;
        g_k0  = Gout.g_k0  ;
        g_co2 = Gout.g_co2 ;
        tmp             = M3d*0         ;
        tmp(iwet(isrf)) = KCO2.*(g_k0.*pco2atm - g_co2)*par.permil ;
        vout.G_dic      = d0(tmp(iwet)) ;

        g_k0_alk  = Gout.g_k0_alk  ;
        g_co2_alk = Gout.g_co2_alk ;
        tmp             = M3d*0         ;
        tmp(iwet(isrf)) = KCO2.*(g_k0_alk.*pco2atm - g_co2_alk)*par.permil ;
        vout.G_alk      = d0(tmp(iwet)) ;

        tmp             = M3d*0         ;
        tmp(iwet(isrf)) = KCO2.*k0*par.permil ;
        vout.G_atm      = tmp(iwet)     ;

        % Hessian
        % [co2surf,k0,Gout] = eqco2(vDICs,vALKs,co2syspar) ;
        % tmp             = M3d*0         ;
        % tmp(iwet(isrf)) = -KCO2.*gg_co2*par.permil ;
        % vout.KGG        = tmp(iwet)     ;    
    end
    %
    if strcmp(Gtype, 'C13')
        %
        pco2atm   = par.pco2atm      ;  % uatm
        vDICs     = par.DIC(isrf)    ;
        vALKs     = par.ALK(isrf)    ;
        co2syspar = par.co2syspar    ;
        scco2 = 2073.1 - 125.62*vSST + 3.6276*vSST.^2 - 0.043219*vSST.^3;
        kw    = kw.*sqrt(660./scco2) ;
        KCO2  = kw/grd.dzt(1)        ;
        %
        % co2surf unit umol/kg;
        % k0 unit mol/kg/atm
        [co2surf,k0] = eqco2(vDICs,vALKs,co2syspar) ;
        co2sat = k0*pco2atm    ; %mol/kg/atm -> umol/kg
        %tmp(iwet(isrf)) = KCO2.*(co2sat - co2surf)*par.permil ;
        %vout.JgDIC = tmp(iwet) ; % umole/kg/s to mmol/m^3/s

        pc13atm  = par.pc13atm     ;    % convert delta c13 to c13 uatm;
        c13surf  = par.DIC13(isrf) ;    % ocean surface c13 concentration  
        c13sat = k0*pc13atm ;           % c13 satuation concentration
      
        tmp    = M3d*0      ;
        alpha_g2dic = par.c13.alpha_g2dic; % gaseous co2 to DIC frationation factor
        alpha_k = par.c13.alpha_k;         % kinectic frationation factor
        alpha_g2aq = par.c13.alpha_g2aq;   % isotopic fractionation factor from gaseous to aqueous CO2
                                           % include fractionation factors in the bulk air-sea flux formula
                                           % see eq (4) in  A. Schmittner et al. (2013) Biogeosciences
        dR13o = par.c13.dR13o(isrf);
        R13o = par.c13.R13o(isrf);
        R13a = par.c13.R13a;
        tmp(iwet(isrf)) = KCO2.*(c13sat.*R13a - (co2surf./alpha_g2dic).*R13o)*alpha_k*alpha_g2aq*par.permil ; 
        vout.JgDIC13 = tmp(iwet) ; % umole/kg/s to mmol/m^3/s
        
        % Gradient
        tmp             = M3d*0         ;
        tmp(iwet(isrf)) = -KCO2.*(co2surf./alpha_g2dic).*dR13o*alpha_k*alpha_g2aq*par.permil;
        vout.G_dic13 = d0(tmp(iwet));        
    end
    %
    if strcmp(Gtype, 'O2')
        KO2  = M3d*0;
        sco2 = 1638.0 - 81.83*vSST + 1.483*vSST.^2 - 0.008004*vSST.^3;
        kw   = kw.*sqrt(660./sco2)      ;
        KO2(iwet(isrf)) = kw/grd.dzt(1) ;
        vout.KO2        = d0(KO2(iwet)) ;
        
        o2sat             = 0*M3d       ;
        o2sat(iwet(isrf)) = 1000*o2sato(vSST,vSSS).*P ; % mmol/m^3
        vout.o2sat        = o2sat(iwet) ;
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
    A0  = 2.00907     ;
    A1  = 3.22014     ;
    A2  = 4.05010     ;
    A3  = 4.94457     ;
    A4  = -2.56847e-1 ;
    A5  =  3.88767    ;
    B0  = -6.24523e-3 ;
    B1  = -7.37614e-3 ;
    B2  = -1.03410e-2 ;
    B3  = -8.17083e-3 ;
    C0  = -4.88682e-7 ;
    TT  = 298.15-T    ;
    TK  = 273.15+T    ;
    TS  = log(TT./TK) ;
    TS2 = TS.^2       ;
    TS3 = TS.^3       ;
    TS4 = TS.^4       ;
    TS5 = TS.^5       ;
    CO  = A0 + A1*TS + A2*TS2 + A3*TS3 + A4*TS4 + A5*TS5+...
          S.*(B0 + B1*TS + B2*TS2 + B3*TS3)+...
          C0*(S.*S)   ;
    o2sat = exp(CO)   ;
    %
    %  Convert from ml/l to mol/m^3
    o2sat = (o2sat/22391.6)*1000.0 ;
end

function [co2,k0,Gout] = eqco2(dic,alk,arg1)
% unpack parameters for co2 system chemistry
    % alk  = arg1.alk   ;
    salt = arg1.salt  ;
    temp = arg1.temp  ;
    pres = arg1.pres  ;
    si   = arg1.si    ;
    po4  = arg1.po4   ;
    
    % co2 system
    a = CO2SYS(abs(alk),abs(dic),1,2,salt,temp,temp,pres,pres,si,po4,1,4,1) ;
    % concentration of co2 in umol/kg
    co2   = a(:,8)    ;
    % co2 solubility k0 
    pco2  = a(:,4)    ; % uatm
    k0    = co2./pco2 ; % mol/kg/atm
    
    % gradient d_dic 
    ax = CO2SYS(abs(real(alk)),abs(real(dic))+sqrt(-1)*eps^3,1,2,salt, ...
                temp,temp,pres,pres,si,po4,1,4,1) ;
    Gout.g_k0  = imag(ax(:,8)./ax(:,4))/eps^3 ;
    Gout.g_co2 = imag(ax(:,8))/eps^3 ;
    clear a ax 
    
    % d_alk
    ax = CO2SYS(abs(real(alk))+sqrt(-1)*eps^3,abs(real(dic)),1,2,salt, ...
                temp,temp,pres,pres,si,po4,1,4,1) ;
    Gout.g_k0_alk  = imag(ax(:,8)./ax(:,4))/eps^3 ;
    Gout.g_co2_alk = imag(ax(:,8))/eps^3 ;
    
    if (nargout>3)
        a = CO2SYS(alk,abs(real(dic))+sqrt(-1)*eps^3,1,2,salt,temp, ...
                   temp,pres,pres,si,po4,1,4,1);
        R      = a(:,14)     ;
        tmp    = R.*co2./dic ;
        gg_co2 = imag(tmp)./eps^3 ;
    end
    % pH (pH units)
    % pH = a(:,3);
end
