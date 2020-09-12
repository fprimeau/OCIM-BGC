%%%%%%%%%%%%%%%%%%%%%%%%  PME part  %%%%%%%%%%%%%%%%%%%%%%%%
function [modT,modS] = PME(par)
M3d = par.M3d;
grd = par.grd;
iwet = par.iwet;
nwet = par.nwet;
TRdiv = par.TRdiv;
dVt = par.dVt;
I    = speye(nwet);   % make an identity matrix;

tmp  = M3d;
tmp(:,:,2:end) = 0;
isrf = find(tmp(iwet)); 
Isrf = d0(tmp(iwet));

B  = TRdiv+Isrf*par.tau_TA;
%
fprintf('Factoring the big matrix for surface salinity 1...'); tic
FB = mfactor(B);
toc;

%
modS = M3d+nan;
Ssurf = par.Salt(iwet);
fprintf('Factoring the big matrix for surface salinity 2...'); tic
modS(iwet) = mfactor(FB,Isrf*Ssurf*par.tau_TA);
toc

%
modT = M3d+nan;
Tsurf = par.Temp(iwet);
fprintf('Factoring the big matrix for surface temperature 1...'); tic
modT(iwet) = mfactor(FB,Isrf*Tsurf*par.tau_TA);
toc
fprintf('\n')
clear memory
clear B 
%
% msk = M3d;
% msk(:,:,2:end) = 0;
% sg   = sum(msk(iwet).*dVt(iwet).*ss(iwet))/sum(msk(iwet).*dVt(iwet));
% test = (ss(iwet)-par.Salt(iwet))/sg;
% pme  = msk(iwet).*test.*(1*par.tau_TA);
% for k = 1:1
    % pme = pme-msk(iwet).*sum(dVt(iwet).*pme)/sum(msk(iwet).*dVt(iwet));
% end

%%%%%%%%%%%%%%%%%%%%%%  END PME part  %%%%%%%%%%%%%%%%%%%%%%%
