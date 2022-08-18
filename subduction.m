addpath('~/Dropbox/myfunc/')
addpath('~/OneDrive - xmu.edu.cn/DATA/'    )
addpath('~/OneDrive - xmu.edu.cn/DATA/OCIM')
load OCIM2_CTL_He.mat
load tempobs_91x180x24.mat
spa = 365.25*24*60*60  ;
% load('~/OneDrive/rDOC/MSK91/CTL_He_PCO_Gamma1to3_POC2DIC_GM15_MODIS_CbPM_aveTeu_diffSig_O2C_uniEta_DOC1e+00_DOP0e+00_xhat.mat')
% load('~/OneDrive/rDOC/MSK91/CTL_He_PCO_Gamma1to3_POC2DIC_GM15_VGPM_aveTeu_diffSig_O2C_uniEta_DOC1e+00_DOP0e+00_xhat.mat')
load('~/OneDrive/rDOC/MSK91/CTL_He_PCO_Gamma1to3_POC2DIC_GM15_CbPM_aveTeu_diffSig_O2C_uniEta_DOC1e+00_DOP0e+00_xhat.mat')
tempobs(tempobs(:)<-2.0) = -2.0 ;
Temp  = tempobs        ;
TRdiv = -output.TR/spa ; %(1/s)
grd   = output.grid    ;
M3d   = output.M3d     ;
iwet  = find(M3d(:))   ;
nwet  = length(iwet)   ;
I     = speye(nwet)    ;
dVt   = grd.DXT3d.*grd.DYT3d.*grd.DZT3d ;

for ji = 1:24
    t2d = Temp(:,:,ji) ;  
    Temp(:,:,ji) = smoothit(grd,M3d,t2d,3,1e5) ;
end 
vT = Temp(iwet) ;
Tz = (vT - min(vT) + 1)./(max(vT) - min(vT)) * 1e-8 ;

Q10C = xhat.Q10C ;
kdC  = xhat.kdC  ;

tf  = (vT - 30)/10 ;
kC  = kdC * Q10C .^ tf ;

tau1 = 1 ./ kC ; % s
L    = M3d ;
L1   = d0(L(iwet) ./ tau1) ;

tau2 = spa/365.25 ;
L(:,:,2:end) = 0 ;
L2 = d0(L(iwet)) / tau2 ;

Z = 0*I ;
M = [[TRdiv+L1,        Z]  ; ...
     [     -L1, TRdiv+L2]] ;

msk = M3d ;
msk(:,:,2:end) = 0  ;
isrf = find(msk(:)) ;
[iy,ix,iz] = find(msk) ;

fprintf('Making the RHS...') ;
tic
iwet = find(M3d) ;
RHS  = zeros(2*length(iwet),length(iy)) ;
for i = 1:length(iy)
    Q = M3d.*0 ;
    Q(iy(i),ix(i),1:2) = 1 ;
    dV = sum(M3d(:).*Q(:).*dVt(:))  ;
    RHS(:,i) = [Q(iwet)./dV;0*iwet] ;
end
toc
fprintf('Factoring the BFM...') ;
tic
FM = mfactor(M') ;
toc
V = [dVt(iwet);dVt(iwet)] ;

fprintf('Solving...') ;
tic
%}
tau = RHS'*mfactor(FM,V)/spa ;
%for i = 1:length(iy)
%  G = mfactor(FM,RHS(:,i));
%  tau(i) = V'*G/(spa);
%  fprintf('%i %f ',i,tau(i)); toc
%end
m = squeeze(M3d(:,:,1))+nan;
for i = 1:length(iy)
  m(iy(i),ix(i)) = tau(i);
end
save mkFigs/doccrestime_SeaWiFS_CbPM.mat m tau