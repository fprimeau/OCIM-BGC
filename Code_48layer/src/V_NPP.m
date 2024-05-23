function [x, a, K, C] = V_NPP(par)
  % vertical distributions of NPP(v_NPP) using three parameters(a, K and C) based on logistic function
  % a: Curver's maximum values, K: steepness of v_NPP, C: sigmoid midpoint of v_NPP

  M3d = par.M3d  ;
  grd = par.grd ;
  nl  = par.nl ;
  ZT_eu = squeeze(grd.ZT3d(1,1,:)) ; 
  ZT_eu(nl+1:end) = [] ;  % It (9) can be changed depending on euphotic zone depth.

  % vertical distributions of NPP(v_NPP) using two parameters(K and C) based on sigmoid function
  % related to integrate NPP from satelite algorithms
  % K: steepness of v_NPP, C: sigmoid midpoint of v_NPP
  K = -0.1 ;
  C =  ZT_eu(4)    ;
  tmp = sum(1 ./ (1+exp(-K*(ZT_eu - C)))) ;
  a = 1/tmp ;
  x = a ./ 1 ./ (1+exp(-K*(ZT_eu - C))) ;
  %
end
