function f = delta1314()
  % addpath ~/MatUtils/seawaterV3.3.1/;
  f.d13 = @(dic13,dic) (dic13./dic/1.12372*1e2-1)*1e3;
  f.d14 = @(dic14,dic,delta13) (dic14./dic/(1.176*1e-12).*(0.975/(1+delta13*1e-3)).^2-1)*1e3;
  f.d14a = @(dic14,dic) -log(dic14./dic)/log(2)*5730;
  f.vm =  @(v,dV,dims) squeeze(nansum(v.*dV,dims))./squeeze(nansum(v.*dV./v,dims));
  f.vs =  @(v,dV,dims) squeeze(nansum(v.*dV,dims))*12/1e18;
  f.nearest = @(x,x0)  find(min(abs(x - x0)) == abs(x - x0));
  R13pdb = 1.12372*1e-2;
  R14oxa = 1.176*1e-12;
  f.d2r13 = @(x) (1+x*1e-3)*R13pdb; % delta c13 to c13/c
  f.D2r14 = @(x,d13) (1+x*1e-3)./((0.975./(1+d13*1e-3)).^2)*R14oxa;
end
