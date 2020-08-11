%%%%%%%%%%%%%%%%%%%%%%%%  PME part  %%%%%%%%%%%%%%%%%%%%%%%%
function [Tmod,Smod] = PME(parm)
M3d = parm.M3d;
grd = parm.grd;
iwet = parm.iwet;
nwet = parm.nwet;
TRdiv = parm.TRdiv;
dVt = parm.dVt;
I    = speye(nwet);   % make an identity matrix;

tmp  = M3d;
tmp(:,:,2:end) = 0;
isrf = find(tmp(iwet)); 
Isrf = d0(tmp(iwet));

B  = TRdiv+Isrf*parm.tau_TA;
%
%fprintf('Factoring the big matrix for surface salinity 1...'); tic
FB = mfactor(B);

%
Smod = M3d+nan;
Ssurf = parm.Salt(iwet);
%fprintf('Factoring the big matrix for surface salinity 2...'); tic
Smod(iwet) = mfactor(FB,Isrf*Ssurf*parm.tau_TA);

%
Tmod = M3d+nan;
Tsurf = parm.Temp(iwet);
%fprintf('Factoring the big matrix for surface temperature 1...'); tic
Tmod(iwet) = mfactor(FB,Isrf*Tsurf*parm.tau_TA);
clear B 
%
% msk = M3d;
% msk(:,:,2:end) = 0;
% sg   = sum(msk(iwet).*dVt(iwet).*ss(iwet))/sum(msk(iwet).*dVt(iwet));
% test = (ss(iwet)-parm.Salt(iwet))/sg;
% pme  = msk(iwet).*test.*(1*parm.tau_TA);
% for k = 1:1
    % pme = pme-msk(iwet).*sum(dVt(iwet).*pme)/sum(msk(iwet).*dVt(iwet));
% end

%%%%%%%%%%%%%%%%%%%%%%  END PME part  %%%%%%%%%%%%%%%%%%%%%%%

end

