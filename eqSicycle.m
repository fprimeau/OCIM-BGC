function [par,Si,Six,Sixx] = eqSicycle(par, x)
on = true; off = false;
M3d =  par.M3d;
iwet = par.iwet;
nwet = par.nwet;
I    = speye(nwet); % make an identity matrix;
par.nx = length(x);
nsx = 0; % count number of Si tunable parameters
% unpack the parameters to be optimized
if (par.opt_alpha == on)
    lalpha = x(par.pindx.lalpha);
    par.alpha  = exp(lalpha);
else
    par.alpha  = par.alpha;
end
% beta
if (par.opt_beta == on)
    lbeta = x(par.pindx.lbeta);
    par.beta  = exp(lbeta);
else
    par.beta  = par.beta;
end
% bsi
if (par.opt_bsi == on)
    nsx = nsx + 1;
    lbsi = x(par.pindx.lbsi);
    par.bsi = exp(lbsi);
else
    par.bsi  = par.bsi;
end
% at
if (par.opt_at == on)
    nsx = nsx + 1;
    lat = x(par.pindx.lat);
    par.at = exp(lat);
else
    par.at  = par.at;
end
% bt
if (par.opt_bt == on)
    nsx = nsx + 1;
    lbt = x(par.pindx.lbt);
    par.bt  = exp(lbt);
else
    par.bt  = par.bt;
end
% aa
if (par.opt_aa == on)
    nsx = nsx + 1;
    par.aa = x(par.pindx.aa);
else
    par.aa  = par.aa;
end
% bb
if (par.opt_bb == on)
    nsx = nsx + 1;
    lbb = x(par.pindx.lbb);
    par.bb  = exp(lbb);
else
    par.bb  = par.bb;
end
par.nsx = nsx;
%
% ++++++++++++++++++++++++++++++++++++++++++++++++++
at  = par.at;
bt  = par.bt;
aa  = par.aa;
bb  = par.bb;
bsi = par.bsi;
kappa_gs = par.kappa_gs;
T   = par.sst(iwet) + 273.15;
kappa_si = at * exp(-bt./T); 
% ++++++++++++++++++++++++++++++++++++++++++++++++++
%
% fixed parameters
SI4 = par.SIL(iwet);
SILbar  = par.SILbar*M3d(iwet);
%
%%%%%%%%%%%%%%%%%%% create Si2P %%%%%%%%%%%%%%%%%%%%%
smsk = M3d;
smsk(:,:,3:end) = 0;
isrf = find(smsk(iwet));
dVs = par.dVt(iwet(isrf));
surface_mean = @(x) sum(x(isrf).*dVs)/sum(dVs);

% compute the mean of the regressor variable
Z = SI4;
mu = surface_mean(Z);
Delta = sqrt(surface_mean((Z-mu).^2));

% standardize the regressor variables
ZR = 0.5+0.5*tanh((Z-mu)/Delta);
% ZR = (Z-min(Z(isrf)))./(max(Z(isrf))-min(Z(isrf)));
%
Si2C = (aa*ZR + bb)./Z;   par.Si2C = Si2C;
dSi2Cdaa = ZR./Z;         par.dSi2Cdaa = dSi2Cdaa;
dSi2Cdbb = 1./Z;          par.dSi2Cdbb = dSi2Cdbb;
%%%%%%%%%%%%%%%%%%%%%%% end %%%%%%%%%%%%%%%%%%%%%%%%
%
PFdiv  = buildPFD(M3d,grd,par,'bSi');
TRdiv = par.TRdiv;
G     = uptake_Si(par);
% build Jacobian matrix
% +++++++++++++++++++++++++++++++++++++++++
% F = ([TRdiv*SIL + G*Si2C*SIL - kappa_si*DSI + kappa_gs*(SIL-SILbar);...
% PFdiv*DSI - G*Si2C*SIL + kappa_si*DSI]);
Jac = [TRdiv+kappa_gs*I+d0(G*Si2C), -d0(kappa_si);...
       -d0(G*Si2C),  PFdiv+d0(kappa_si)];

RHS = [SILbar*kappa_gs; ...
       sparse(nwet,1)];

FD = mfactor(Jac);
Si = mfactor(FD,RHS);
SIL = Si(1:nwet);
DSI = Si(nwet+1:end);
%% +++++++++++++++++++++++++++++++++++++++++++
if nargout > 1
    [~,Gx] = uptake_Si(par);
    % sigma
    if (par.opt_sigma == on)
        tmp =  [-d0(Gx(:,par.pindx.lsigma))*SIL.*Si2C;...
                d0(Gx(:,par.pindx.lsigma))*SIL.*Si2C]   ;

        Six(:,par.pindx.lsigma) = mfactor(FD, tmp);
    end
    
    % interpp
    if (par.opt_interpp == on)
        tmp = [-d0(Gx(:,par.pindx.linterpp))*SIL.*Si2C ;...
               d0(Gx(:,par.pindx.linterpp))*SIL.*Si2C]   ;

        Six(:,par.pindx.linterpp) = mfactor(FD, tmp);
    end
    
    % slopep
    if (par.opt_slopep == on)
        tmp = [-d0(Gx(:,par.pindx.slopep))*SIL.*Si2C ;...
               d0(Gx(:,par.pindx.slopep))*SIL.*Si2C];

        Six(:,par.pindx.slopep) = mfactor(FD, tmp);
    end
    
    % kappa_dp
    if (par.opt_kappa_dp == on)
        tmp = [-d0(Gx(:,par.pindx.lkappa_dp))*SIL.*Si2C;...
               d0(Gx(:,par.pindx.lkappa_dp))*SIL.*Si2C];

        Six(:,par.pindx.lkappa_dp) = mfactor(FD, tmp);
    end

    % alpha
    if (par.opt_alpha == on)
        tmp = [-d0(Gx(:,par.pindx.lalpha))*SIL.*Si2C; ...
               d0(Gx(:,par.pindx.lalpha))*SIL.*Si2C];

        Six(:,par.pindx.lalpha) = mfactor(FD, tmp);
    end

    % beta
    if (par.opt_beta == on)
        tmp = [-d0(Gx(:,par.pindx.lbeta))*SIL.*Si2C; ...
               d0(Gx(:,par.pindx.lbeta))*SIL.*Si2C];

        Six(:,par.pindx.lbeta) = mfactor(FD, tmp);
    end

    % bsi
    if (par.opt_bsi == on)
        [~,gout] = buildPFD(M3d,grd,par,'bSi');
        dPFDdb = gout.PFD_b;
        tmp = bsi*[DSI*0;  -dPFDdb*DSI];

        Six(:,par.pindx.lbsi) = mfactor(FD, tmp);
    end
    
    % at
    if (par.opt_at == on)
        [~,gout] = buildPFD(M3d,grd,par,'bSi');
        dPFDdat = gout.PFD_at;
        dkdat = exp(-bt./T);
        tmp = at*[dkdat.*DSI; ...
                  -(d0(dkdat)+dPFDdat)*DSI];

        Six(:,par.pindx.lat) = mfactor(FD, tmp);
    end

    % bt
    if (par.opt_bt == on)
        [~,gout] = buildPFD(M3d,grd,par,'bSi');
        dPFDdbt = gout.PFD_bt;
        dkdbt = -(at*exp(-bt./T))./T;
        tmp = bt*[dkdbt.*DSI; ...
                  -(d0(dkdbt)+dPFDdbt)*DSI];

        Six(:,par.pindx.lbt) = mfactor(FD, tmp);
    end

    % aa
    if (par.opt_aa == on)
        tmp = [-d0(G*SIL)*dSi2Cdaa; ...
               d0(G*SIL)*dSi2Cdaa];

        Six(:,par.pindx.aa) = mfactor(FD, tmp);
    end

    % bb
    if (par.opt_bb == on)
        tmp = bb*[-d0(G*SIL)*dSi2Cdbb; ...
                  d0(G*SIL)*dSi2Cdbb];

        Six(:,par.pindx.lbb) = mfactor(FD, tmp);
    end

end
SILx = Six(1:nwet,:);
DSIx = Six(nwet+1:end,:);
%% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
if (nargout > 2)
    nx = par.npx + par.nsx;
    ncs = nchoosek(nx,2)+nx;
    Sixx = sparse(2*nwet,ncs);
    % Compute the hessian of the solution wrt the parameters
    DIP = par.DIP;
    DIPx  = par.Px(1:nwet,:);
    DIPxx = par.Pxx(1:nwet,:);
    [~,~,Gxx,Gp,par] = uptake_Si(par);
    dGpdalpha = par.dGpdalpha;
    dGpdbeta  = par.dGpdbeta;
    % sigma sigma
    kk = 1;
    if (par.opt_sigma)
        tmp = [-d0(Gxx(:,kk))*SIL.*Si2C; ...
               d0(Gxx(:,kk))*SIL.*Si2C] + ...
              2*[-d0(Gx(:,par.pindx.lsigma))*SILx(:,par.pindx.lsigma).*Si2C;...
                 d0(Gx(:,par.pindx.lsigma))*SILx(:,par.pindx.lsigma).*Si2C];

        Sixx(:,kk) = mfactor(FD, tmp);
        kk = kk + 1;
    end
    % sigma kappa_dp
    if (par.opt_sigma & par.opt_kappa_dp)
        tmp = [-d0(Gxx(:,kk))*SIL.*Si2C; ...
               d0(Gxx(:,kk))*SIL.*Si2C] + ...
              [-d0(Gx(:,par.pindx.lsigma))*SILx(:,par.pindx.lkappa_dp).*Si2C;...
               d0(Gx(:,par.pindx.lsigma))*SILx(:,par.pindx.lkappa_dp).*Si2C] + ...
              [-d0(Gx(:,par.pindx.lkappa_dp))*SILx(:,par.pindx.lsigma).*Si2C;...
               d0(Gx(:,par.pindx.lkappa_dp))*SILx(:,par.pindx.lsigma).*Si2C];

        Sixx(:,kk) = mfactor(FD, tmp);
        kk = kk + 1;
    end
    % sigma slopep
    if (par.opt_sigma & par.opt_slopep)
        tmp = [-d0(Gxx(:,kk))*SIL.*Si2C; ...
               d0(Gxx(:,kk))*SIL.*Si2C] + ...
              [-d0(Gx(:,par.pindx.lsigma))*SILx(:,par.pindx.slopep).*Si2C;...
               d0(Gx(:,par.pindx.lsigma))*SILx(:,par.pindx.slopep).*Si2C] + ...
              [-d0(Gx(:,par.pindx.slopep))*SILx(:,par.pindx.lsigma).*Si2C;...
               d0(Gx(:,par.pindx.slopep))*SILx(:,par.pindx.lsigma).*Si2C];

        Sixx(:,kk) = mfactor(FD, tmp);
        kk = kk + 1;
    end
    % sigma interpp
    if (par.opt_sigma == on & par.opt_interpp == on)
        tmp = [-d0(Gxx(:,kk))*SIL.*Si2C; ...
               d0(Gxx(:,kk))*SIL.*Si2C] + ...
              [-d0(Gx(:,par.pindx.lsigma))*SILx(:,par.pindx.linterpp).*Si2C;...
               d0(Gx(:,par.pindx.lsigma))*SILx(:,par.pindx.linterpp).*Si2C]+...
              [-d0(Gx(:,par.pindx.linterpp))*SILx(:,par.pindx.lsigma).*Si2C;...
               d0(Gx(:,par.pindx.linterpp))*SILx(:,par.pindx.lsigma).*Si2C];

        Sixx(:,kk) = mfactor(FD, tmp);
        kk = kk + 1;
    end
    % sigma alpha
    if (par.opt_sigma == on & par.opt_alpha == on)
        tmp = [-d0(Gxx(:,kk))*SIL.*Si2C; ...
               d0(Gxx(:,kk))*SIL.*Si2C] + ...
              [-d0(Gx(:,par.pindx.lsigma))*SILx(:,par.pindx.lalpha).*Si2C;...
               d0(Gx(:,par.pindx.lsigma))*SILx(:,par.pindx.lalpha).*Si2C]+...
              [-d0(Gx(:,par.pindx.lalpha))*SILx(:,par.pindx.lsigma).*Si2C;...
               d0(Gx(:,par.pindx.lalpha))*SILx(:,par.pindx.lsigma).*Si2C];
        
        Sixx(:,kk) = mfactor(FD, tmp);
        kk = kk + 1;
    end
    % sigma beta
    if (par.opt_sigma == on & par.opt_beta == on)
        tmp = [-d0(Gxx(:,kk))*SIL.*Si2C; ...
               d0(Gxx(:,kk))*SIL.*Si2C] + ...
              [-d0(Gx(:,par.pindx.lsigma))*SILx(:,par.pindx.lbeta).*Si2C;...
               d0(Gx(:,par.pindx.lsigma))*SILx(:,par.pindx.lbeta).*Si2C]+...
              [-d0(Gx(:,par.pindx.lbeta))*SILx(:,par.pindx.lsigma).*Si2C;...
               d0(Gx(:,par.pindx.lbeta))*SILx(:,par.pindx.lsigma).*Si2C];
        
        Sixx(:,kk) = mfactor(FD, tmp);
        kk = kk + 1;
    end
    % kappa_dp kappa_dp
    if (par.opt_kappa_dp == on)
        tmp = [-d0(Gxx(:,kk))*SIL.*Si2C; ...
               d0(Gxx(:,kk))*SIL.*Si2C] + ...
              2*[-d0(Gx(:,par.pindx.lkappa_dp))*SILx(:,par.pindx.lkappa_dp).*Si2C;...
                 d0(Gx(:,par.pindx.lkappa_dp))*SILx(:,par.pindx.lkappa_dp).*Si2C];

        Sixx(:,kk) = mfactor(FD, tmp);
        kk = kk + 1;
    end
    % kappa_dp slopep
    if (par.opt_kappa_dp == on & par.opt_slopep == on)
        tmp = [-d0(Gxx(:,kk))*SIL.*Si2C; ...
               d0(Gxx(:,kk))*SIL.*Si2C] + ...
              [-d0(Gx(:,par.pindx.lkappa_dp))*SILx(:,par.pindx.slopep).*Si2C;...
               d0(Gx(:,par.pindx.lkappa_dp))*SILx(:,par.pindx.slopep).*Si2C]+...
              [-d0(Gx(:,par.pindx.slopep))*SILx(:,par.pindx.lkappa_dp).*Si2C;...
               d0(Gx(:,par.pindx.slopep))*SILx(:,par.pindx.lkappa_dp).*Si2C];

        Sixx(:,kk) = mfactor(FD, tmp);
        kk = kk + 1;
    end
    % kappa_dp interpp
    if (par.opt_kappa_dp == on & par.opt_interpp == on)
        tmp = [-d0(Gxx(:,kk))*SIL.*Si2C; ...
               d0(Gxx(:,kk))*SIL.*Si2C] + ...
              [-d0(Gx(:,par.pindx.lkappa_dp))*SILx(:,par.pindx.linterpp).*Si2C;...
               d0(Gx(:,par.pindx.lkappa_dp))*SILx(:,par.pindx.linterpp).*Si2C]+...
              [-d0(Gx(:,par.pindx.linterpp))*SILx(:,par.pindx.lkappa_dp).*Si2C;...
               d0(Gx(:,par.pindx.linterpp))*SILx(:,par.pindx.lkappa_dp).*Si2C];

        Sixx(:,kk) = mfactor(FD, tmp);
        kk = kk + 1;
    end
    % kappa_dp alpha
    if (par.opt_kappa_dp == on & par.opt_alpha == on)
        tmp = [-d0(Gxx(:,kk))*SIL.*Si2C; ...
               d0(Gxx(:,kk))*SIL.*Si2C] + ...
              [-d0(Gx(:,par.pindx.lkappa_dp))*SILx(:,par.pindx.lalpha).*Si2C;...
               d0(Gx(:,par.pindx.lkappa_dp))*SILx(:,par.pindx.lalpha).*Si2C]+...
              [-d0(Gx(:,par.pindx.lalpha))*SILx(:,par.pindx.lkappa_dp).*Si2C;...
               d0(Gx(:,par.pindx.lalpha))*SILx(:,par.pindx.lkappa_dp).*Si2C];
        
        Sixx(:,kk) = mfactor(FD, tmp);
        kk = kk + 1;
    end
    % kappa_dp beta
    if (par.opt_kappa_dp == on & par.opt_beta == on)
        tmp = [-d0(Gxx(:,kk))*SIL.*Si2C; ...
               d0(Gxx(:,kk))*SIL.*Si2C] + ...
              [-d0(Gx(:,par.pindx.lkappa_dp))*SILx(:,par.pindx.lbeta).*Si2C;...
               d0(Gx(:,par.pindx.lkappa_dp))*SILx(:,par.pindx.lbeta).*Si2C]+...
              [-d0(Gx(:,par.pindx.lbeta))*SILx(:,par.pindx.lkappa_dp).*Si2C;...
               d0(Gx(:,par.pindx.lbeta))*SILx(:,par.pindx.lkappa_dp).*Si2C];
        
        Sixx(:,kk) = mfactor(FD, tmp);
        kk = kk + 1;
    end
    % slopep slopep
    if (par.opt_slopep == on)
        tmp = [-d0(Gxx(:,kk))*SIL.*Si2C; ...
               d0(Gxx(:,kk))*SIL.*Si2C] + ...
              2*[-d0(Gx(:,par.pindx.slopep))*SILx(:,par.pindx.slopep).*Si2C;...
                 d0(Gx(:,par.pindx.slopep))*SILx(:,par.pindx.slopep).*Si2C];
        
        Sixx(:,kk) = mfactor(FD, tmp);
        kk = kk + 1;
    end
    % slopep interpp
    if (par.opt_slopep == on & par.opt_interpp == on)
        tmp = [-d0(Gxx(:,kk))*SIL.*Si2C; ...
               d0(Gxx(:,kk))*SIL.*Si2C] + ...
              [-d0(Gx(:,par.pindx.slopep))*SILx(:,par.pindx.linterpp).*Si2C;...
               d0(Gx(:,par.pindx.slopep))*SILx(:,par.pindx.linterpp).*Si2C]+...
              [-d0(Gx(:,par.pindx.linterpp))*SILx(:,par.pindx.slopep).*Si2C;...
               d0(Gx(:,par.pindx.linterpp))*SILx(:,par.pindx.slopep).*Si2C];

        Sixx(:,kk) = mfactor(FD, tmp);
        kk = kk + 1;
    end
    % slopep alpha
    if (par.opt_slopep == on & par.opt_alpha == on)
        tmp = [-d0(Gxx(:,kk))*SIL.*Si2C; ...
               d0(Gxx(:,kk))*SIL.*Si2C] + ...
              [-d0(Gx(:,par.pindx.slopep))*SILx(:,par.pindx.lalpha).*Si2C;...
               d0(Gx(:,par.pindx.slopep))*SILx(:,par.pindx.lalpha).*Si2C]+...
              [-d0(Gx(:,par.pindx.lalpha))*SILx(:,par.pindx.slopep).*Si2C;...
               d0(Gx(:,par.pindx.lalpha))*SILx(:,par.pindx.slopep).*Si2C];

        Sixx(:,kk) = mfactor(FD, tmp);
        kk = kk + 1;
    end
    % slopep beta
    if (par.opt_slopep == on & par.opt_beta == on)
        tmp = [-d0(Gxx(:,kk))*SIL.*Si2C; ...
               d0(Gxx(:,kk))*SIL.*Si2C] + ...
              [-d0(Gx(:,par.pindx.slopep))*SILx(:,par.pindx.lbeta).*Si2C;...
               d0(Gx(:,par.pindx.slopep))*SILx(:,par.pindx.lbeta).*Si2C]+...
              [-d0(Gx(:,par.pindx.lbeta))*SILx(:,par.pindx.slopep).*Si2C;...
               d0(Gx(:,par.pindx.lbeta))*SILx(:,par.pindx.slopep).*Si2C];

        Sixx(:,kk) = mfactor(FD, tmp);
        kk = kk + 1;
    end
    % interpp interpp
    if (par.opt_interpp == on)
        tmp = [-d0(Gxx(:,kk))*SIL.*Si2C; ...
               d0(Gxx(:,kk))*SIL.*Si2C] + ...
              2*[-d0(Gx(:,par.pindx.linterpp))*SILx(:,par.pindx.linterpp).*Si2C;...
                 d0(Gx(:,par.pindx.linterpp))*SILx(:,par.pindx.linterpp).*Si2C];

        Sixx(:,kk) = mfactor(FD, tmp);
        kk = kk + 1;
    end
    % interpp alpha
    if (par.opt_interpp == on & par.opt_alpha == on)
        tmp = [-d0(Gxx(:,kk))*SIL.*Si2C; ...
               d0(Gxx(:,kk))*SIL.*Si2C] + ...
              [-d0(Gx(:,par.pindx.linterpp))*SILx(:,par.pindx.lalpha).*Si2C; ...
               d0(Gx(:,par.pindx.linterpp))*SILx(:,par.pindx.lalpha).*Si2C]+ ...
              [-d0(Gx(:,par.pindx.lalpha))*SILx(:,par.pindx.linterpp).*Si2C; ...
               d0(Gx(:,par.pindx.lalpha))*SILx(:,par.pindx.linterpp).*Si2C];

        Sixx(:,kk) = mfactor(FD, tmp);
        kk = kk + 1;
    end
    % interpp beta
    if (par.opt_interpp == on & par.opt_beta == on)
        tmp = [-d0(Gxx(:,kk))*SIL.*Si2C; ...
               d0(Gxx(:,kk))*SIL.*Si2C] + ...
              [-d0(Gx(:,par.pindx.linterpp))*SILx(:,par.pindx.lbeta).*Si2C; ...
               d0(Gx(:,par.pindx.linterpp))*SILx(:,par.pindx.lbeta).*Si2C]+ ...
              [-d0(Gx(:,par.pindx.lbeta))*SILx(:,par.pindx.linterpp).*Si2C; ...
               d0(Gx(:,par.pindx.lbeta))*SILx(:,par.pindx.linterpp).*Si2C];

        Sixx(:,kk) = mfactor(FD, tmp);
        kk = kk + 1;
    end
    % alpha alpha
    if (par.opt_alpha == on)
        tmp = [-d0(Gxx(:,kk))*SIL.*Si2C; ...
               d0(Gxx(:,kk))*SIL.*Si2C] + ...
              2*[-d0(Gx(:,par.pindx.lalpha))*SILx(:,par.pindx.lalpha).*Si2C; ...
                 d0(Gx(:,par.pindx.lalpha))*SILx(:,par.pindx.lalpha).*Si2C];

        Sixx(:,kk) = mfactor(FD, tmp);
        kk = kk + 1;
    end
    % alpha beta
    if (par.opt_alpha == on & par.opt_beta == on)
        tmp = [-d0(Gxx(:,kk))*SIL.*Si2C; ...
               d0(Gxx(:,kk))*SIL.*Si2C] + ...
              [-d0(Gx(:,par.pindx.lalpha))*SILx(:,par.pindx.lbeta).*Si2C; ...
               d0(Gx(:,par.pindx.lalpha))*SILx(:,par.pindx.lbeta).*Si2C]+ ...
              [-d0(Gx(:,par.pindx.lbeta))*SILx(:,par.pindx.lalpha).*Si2C; ...
               d0(Gx(:,par.pindx.lbeta))*SILx(:,par.pindx.lalpha).*Si2C];

        Sixx(:,kk) = mfactor(FD, tmp);
        kk = kk + 1;
    end
    % beta beta
    if (par.opt_beta == on)
        tmp = [-d0(Gxx(:,kk))*SIL.*Si2C; ...
               d0(Gxx(:,kk))*SIL.*Si2C] + ...
              2*[-d0(Gx(:,par.pindx.lbeta))*SILx(:,par.pindx.lbeta).*Si2C; ...
                 d0(Gx(:,par.pindx.lbeta))*SILx(:,par.pindx.lbeta).*Si2C];

        Sixx(:,kk) = mfactor(FD, tmp);
        kk = kk + 1;
    end
    % --------------------------------------------------------
    % Simodel parameters
    % sigma bsi
    if (par.opt_sigma == on & par.opt_bsi == on)
        tmp = bsi*[0*DSI; ...
                   -dPFDdb*DSIx(:,par.pindx.lsigma)] + ...
              [-d0(Gx(:,par.pindx.lsigma))*SILx(:,par.pindx.lbsi).*Si2C; ...
               d0(Gx(:,par.pindx.lsigma))*SILx(:,par.pindx.lbsi).*Si2C];
              
        Sixx(:,kk) = mfactor(FD, tmp);
        kk = kk + 1;
    end
    % sigma at
    if (par.opt_sigma == on & par.opt_at == on)
        tmp = at*[d0(dkdat)*DSIx(:,par.pindx.lsigma); ...
                  -(d0(dkdat) + dPFDdat)*DSIx(:,par.pindx.lsigma)]+ ...
              [-d0(Gx(:,par.pindx.lsigma))*SILx(:,par.pindx.lat).*Si2C; ...
               d0(Gx(:,par.pindx.lsigma))*SILx(:,par.pindx.lat).*Si2C];

        Sixx(:,kk) = mfactor(FD, tmp);
        kk = kk + 1;
    end
    % sigma bt
    if (par.opt_sigma == on & par.opt_bt == on)
        tmp = bt*[d0(dkdbt)*DSIx(:,par.pindx.lsigma); ...
                  -(d0(dkdbt)+dPFDdbt)*DSIx(:,par.pindx.lsigma)]+ ...
              [-d0(Gx(:,par.pindx.lsigma))*SILx(:,par.pindx.lbt).*Si2C; ...
               d0(Gx(:,par.pindx.lsigma))*SILx(:,par.pindx.lbt).*Si2C];

        Sixx(:,kk) = mfactor(FD, tmp);
        kk = kk + 1;
    end
    % sigma aa
    if (par.opt_sigma == on & par.opt_aa == on)
        tmp = [-d0(Gxx(:,kk))*SIL; ...
               d0(Gxx(:,kk))*SIL] + ...
              [-d0(G*SILx(:,par.pindx.lsigma))*dSi2Cdaa; ...
               d0(G*SILx(:,par.pindx.lsigma))*dSi2Cdaa]+ ...
              [-d0(Gx(:,par.pindx.lsigma))*SILx(:,par.pindx.aa).*Si2C; ...
               d0(Gx(:,par.pindx.lsigma))*SILx(:,par.pindx.aa).*Si2C];

        Sixx(:,kk) = mfactor(FD, tmp);
        kk = kk + 1;
    end
    % sigma bb
    if (par.opt_sigma == on & par.opt_bb == on)
        tmp = [-d0(Gxx(:,kk))*SIL; ...
               d0(Gxx(:,kk))*SIL] + ...
              bb*[-d0(G*SILx(:,par.pindx.lsigma))*dSi2Cdbb; ...
                  d0(G*SILx(:,par.pindx.lsigma))*dSi2Cdbb]+ ...
              [-d0(Gx(:,par.pindx.lsigma))*SILx(:,par.pindx.lbb).*Si2C; ...
               d0(Gx(:,par.pindx.lsigma))*SILx(:,par.pindx.lbb).*Si2C];

        Sixx(:,kk) = mfactor(FD, tmp);
        kk = kk + 1;
    end
    % kappa_dp bsi
    if (par.opt_kappa_dp == on & par.opt_bsi == on)
        tmp = bsi*[0*DSI; ...
                   -(dPFDdb)*DSIx(:,par.pindx.lkappa_dp)]+ ...
              [-d0(Gx(:,par.pindx.lkappa_dp))*SILx(:,par.pindx.lbsi).*Si2C; ...
               d0(Gx(:,par.pindx.lkappa_dp))*SILx(:,par.pindx.lbsi).*Si2C];

        Sixx(:,kk) = mfactor(FD, tmp);
        kk = kk + 1;
    end
    % kappa_dp at
    if (par.opt_kappa_dp == on & par.opt_at == on)
        tmp = at*[d0(dkdat)*DSIx(:,par.pindx.lkappa_dp); ...
                  -(d0(dkdat)+dPFDdat)*DSIx(:,par.pindx.lkappa_dp)]+ ...
              [-d0(Gx(:,par.pindx.lkappa_dp))*SILx(:,par.pindx.lat).*Si2C; ...
               d0(Gx(:,par.pindx.lkappa_dp))*SILx(:,par.pindx.lat).*Si2C];

        Sixx(:,kk) = mfactor(FD, tmp);
        kk = kk + 1;
    end
    % kappa_dp bt
    if (par.opt_kappa_dp == on & par.opt_bt == on)
        tmp = bt*[d0(dkdbt)*DSIx(:,par.pindx.lkappa_dp); ...
                  -(d0(dkdbt)+dPFDdbt)*DSIx(:,par.pindx.lkappa_dp)]+ ...
              [-d0(Gx(:,par.pindx.lkappa_dp))*SILx(:,par.pindx.lbt).*Si2C; ...
               d0(Gx(:,par.pindx.lkappa_dp))*SILx(:,par.pindx.lbt).*Si2C];

        Sixx(:,kk) = mfactor(FD, tmp); 
        kk = kk + 1;
    end
    % kappa_dp aa
    if (par.opt_kappa_dp == on & par.opt_aa == on)
        tmp =  [-d0(Gxx(:,kk))*SIL; ...
                d0(Gxx(:,kk))*SIL] + ...
               [-d0(G*SILx(:,par.pindx.lkappa_dp))*dSi2Cdaa; ...
                d0(G*SILx(:,par.pindx.lkappa_dp))*dSi2Cdaa]+ ...
               [-d0(Gx(:,par.pindx.lkappa_dp))*SILx(:,par.pindx.aa).*Si2C; ...
                d0(Gx(:,par.pindx.lkappa_dp))*SILx(:,par.pindx.aa).*Si2C];

        Sixx(:,kk) = mfactor(FD, tmp);
        kk = kk + 1;
    end
    % kappa_dp bb
    if (par.opt_kappa_dp == on & par.opt_bb == on)
        tmp = [-d0(Gxx(:,kk))*SIL; ...
               d0(Gxx(:,kk))*SIL] + ...
              bb*[-d0(G*SILx(:,par.pindx.lkappa_dp))*dSi2Cdbb; ...
                  d0(G*SILx(:,par.pindx.lkappa_dp))*dSi2Cdbb]+ ...
              [-d0(Gx(:,par.pindx.lkappa_dp))*SILx(:,par.pindx.lbb).*Si2C; ...
               d0(Gx(:,par.pindx.lkappa_dp))*SILx(:,par.pindx.lbb).*Si2C];

        Sixx(:,kk) = mfactor(FD, tmp);
        kk = kk + 1;
    end
    % slope bsi
    if (par.opt_slopep == on & par.opt_bsi == on)
        tmp = bsi*[0*DSI; ...
                   -dPFDdb*DSIx(:,par.pindx.slopep)]+ ...
              [-d0(Gx(:,par.pindx.slopep))*SILx(:,par.pindx.lbsi).*Si2C; ...
               d0(Gx(:,par.pindx.slopep))*SILx(:,par.pindx.lbsi).*Si2C];

        Sixx(:,kk) = mfactor(FD, tmp);
        kk = kk + 1;
    end
    % slope at
    if (par.opt_slopep == on & par.opt_at == on)
        tmp = at*[d0(dkdat)*DSIx(:,par.pindx.slopep); ...
                  -(d0(dkdat)+dPFDdat)*DSIx(:,par.pindx.slopep)]+ ...
              [-d0(Gx(:,par.pindx.slopep))*SILx(:,par.pindx.lat).*Si2C; ...
               d0(Gx(:,par.pindx.slopep))*SILx(:,par.pindx.lat).*Si2C];

        Sixx(:,kk) = mfactor(FD, tmp);
        kk = kk + 1;
    end
    % slope bt
    if (par.opt_slopep == on & par.opt_bt == on)
        tmp = bt*[d0(dkdbt)*DSIx(:,par.pindx.slopep); ...
                  -(d0(dkdbt)+dPFDdbt)*DSIx(:,par.pindx.slopep)]+ ...
              [-d0(Gx(:,par.pindx.slopep))*SILx(:,par.pindx.lbt).*Si2C; ...
               d0(Gx(:,par.pindx.slopep))*SILx(:,par.pindx.lbt).*Si2C];

        Sixx(:,kk) = mfactor(FD, tmp);
        kk = kk + 1;
    end
    % slopep aa
    if (par.opt_slopep == on & par.opt_aa == on)
        tmp = [-d0(Gxx(:,kk))*SIL; ...
               d0(Gxx(:,kk))*SIL] + ...
              [-d0(G*SILx(:,par.pindx.slopep))*dSi2Cdaa; ...
               d0(G*SILx(:,par.pindx.slopep))*dSi2Cdaa] + ...
              [-d0(Gx(:,par.pindx.slopep))*SILx(:,par.pindx.aa).*Si2C;  ...
               d0(Gx(:,par.pindx.slopep))*SILx(:,par.pindx.aa).*Si2C];

        Sixx(:,kk) = mfactor(FD, tmp);
        kk = kk + 1;
    end
    % slopep bb
    if (par.opt_slopep == on & par.opt_bb == on)
        tmp = [-d0(Gxx(:,kk))*SIL; ...
               d0(Gxx(:,kk))*SIL] + ...
              bb*[-d0(G*SILx(:,par.pindx.slopep))*dSi2Cdbb; ...
                  d0(G*SILx(:,par.pindx.slopep))*dSi2Cdbb] + ...
              [-d0(Gx(:,par.pindx.slopep))*SILx(:,par.pindx.lbb).*Si2C; ...
               d0(Gx(:,par.pindx.slopep))*SILx(:,par.pindx.lbb).*Si2C];

        Sixx(:,kk) = mfactor(FD, tmp);
        kk = kk + 1;
    end
    % interpp bsi
    if (par.opt_interpp == on & par.opt_bsi == on)
        tmp = bsi*[0*DSI; ...
                   -dPFDdb*DSIx(:, par.pindx.linterpp)] + ...
              [-d0(Gx(:,par.pindx.linterpp))*SILx(:,par.pindx.lbsi).*Si2C; ...
               d0(Gx(:,par.pindx.linterpp))*SILx(:,par.pindx.lbsi).*Si2C];

        Sixx(:,kk) = mfactor(FD, tmp);
        kk = kk + 1;
    end
    % interpp at
    if (par.opt_interpp == on & par.opt_at == on)
        tmp = at*[d0(dkdat)*DSIx(:, par.pindx.linterpp); ...
                  -(d0(dkdat)+dPFDdat)*DSIx(:, par.pindx.linterpp)]+ ...
              [-d0(Gx(:,par.pindx.linterpp))*SILx(:,par.pindx.lat).*Si2C; ...
               d0(Gx(:,par.pindx.linterpp))*SILx(:,par.pindx.lat).*Si2C];

        Sixx(:,kk) = mfactor(FD, tmp);
        kk = kk + 1;
    end
    % interpp bt
    if (par.opt_interpp == on & par.opt_bt == on)
        tmp = bt*[d0(dkdbt)*DSIx(:, par.pindx.linterpp); ...
                  -(d0(dkdbt)+dPFDdbt)*DSIx(:,par.pindx.linterpp)]+ ...
              [-d0(Gp*DIPx(:,par.pindx.linterpp))*SILx(:,par.pindx.lbt).*Si2C; ...
               d0(Gp*DIPx(:,par.pindx.linterpp))*SILx(:,par.pindx.lbt).*Si2C];

        Sixx(:,kk) = mfactor(FD, tmp);
        kk = kk + 1;
    end
    % interpp aa
    if (par.opt_interpp == on & par.opt_aa == on)
        tmp = [-d0(Gxx(:,kk))*SIL; ...
               d0(Gxx(:,kk))*SIL] + ...
              [-d0(G*SILx(:,par.pindx.linterpp))*dSi2Cdaa; ...
               d0(G*SILx(:,par.pindx.linterpp))*dSi2Cdaa] + ...
              [-d0(Gx(:,par.pindx.linterpp))*SILx(:,par.pindx.aa).*Si2C; ...
               d0(Gx(:,par.pindx.linterpp))*SILx(:,par.pindx.aa).*Si2C];

        Sixx(:,kk) = mfactor(FD, tmp);
        kk = kk + 1;
    end
    % interpp bb
    if (par.opt_interpp == on & par.opt_bb == on)
        tmp = [-d0(Gxx(:,kk))*SIL; ...
               d0(Gxx(:,kk))*SIL] + ...
              bb*[-d0(G*SILx(:,par.pindx.linterpp))*dSi2Cdbb; ...
                  d0(G*SILx(:,par.pindx.linterpp))*dSi2Cdbb]+ ...
              [-d0(Gx(:,par.pindx.linterpp))*SILx(:,par.pindx.lbb).*Si2C; ...
               d0(Gx(:,par.pindx.linterpp))*SILx(:,par.pindx.lbb).*Si2C];

        Sixx(:,kk) = mfactor(FD, tmp);
        kk = kk + 1;
    end
    % alpha bsi
    if (par.opt_alpha == on & par.opt_bsi == on)
        tmp = bsi*[0*DSI; ...
                   -dPFDdb*DSIx(:,par.pindx.lalpha)] + ...
              [-d0(Gx(:,par.pindx.lalpha))*SILx(:,par.pindx.lbsi).*Si2C; ...
               d0(Gx(:,par.pindx.lalpha))*SILx(:,par.pindx.lbsi).*Si2C];

        Sixx(:,kk) = mfactor(FD, tmp);
        kk = kk + 1;
    end
    % alpha at
    if (par.opt_alpha == on & par.opt_at == on)
        tmp = at*[d0(dkdat)*DSIx(:,par.pindx.lalpha); ...
                  -(d0(dkdat)+dPFDdat)*DSIx(:,par.pindx.lalpha)] + ...
              [-d0(Gx(:,par.pindx.lalpha))*SILx(:,par.pindx.lat).*Si2C; ...
               d0(Gx(:,par.pindx.lalpha))*SILx(:,par.pindx.lat).*Si2C];

        Sixx(:,kk) = mfactor(FD, tmp);
        kk = kk + 1;
    end
    % alpha bt
    if (par.opt_alpha == on & par.opt_bt == on)
        tmp = bt*[d0(dkdbt)*DSIx(:,par.pindx.lalpha); ...
                  -(d0(dkdbt) + dPFDdbt)*DSIx(:,par.pindx.lalpha)] + ...
              [-d0(Gx(:,par.pindx.lalpha))*SILx(:,par.pindx.lbt).*Si2C; ...
               d0(Gx(:,par.pindx.lalpha))*SILx(:,par.pindx.lbt).*Si2C];

        Sixx(:,kk) = mfactor(FD, tmp);
        kk = kk + 1;
    end
    % alpha aa
    if (par.opt_alpha == on & par.opt_aa == on)
        tmp = [-d0(Gxx(:,kk))*SIL; ...
               d0(Gxx(:,kk))*SIL] + ...
              [-d0(G*SILx(:,par.pindx.lalpha))*dSi2Cdaa; ...
               d0(G*SILx(:,par.pindx.lalpha))*dSi2Cdaa] + ...
              [-d0(Gx(:,par.pindx.lalpha))*SILx(:,par.pindx.aa).*Si2C; ...
               d0(Gx(:,par.pindx.lalpha))*SILx(:,par.pindx.aa).*Si2C];

        Sixx(:,kk) = mfactor(FD, tmp);
        kk = kk + 1;
    end
    % alpha bb
    if (par.opt_alpha == on & par.opt_bb == on)
        tmp = [-d0(Gxx(:,kk))*SIL; ...
               d0(Gxx(:,kk))*SIL] + ...
              bb*[-d0(G*SILx(:,par.pindx.lalpha))*dSi2Cdbb; ...
                  d0(G*SILx(:,par.pindx.lalpha))*dSi2Cdbb] + ...
              [-d0(Gx(:,par.pindx.lalpha))*SILx(:,par.pindx.lbb).*Si2C; ...
               d0(Gx(:,par.pindx.lalpha))*SILx(:,par.pindx.lbb).*Si2C];

        Sixx(:,kk) = mfactor(FD, tmp);
        kk = kk + 1;
    end
    % beta bsi
    if (par.opt_beta == on & par.opt_bsi == on)
        tmp = bsi*[0*DSI; ...
                   -dPFDdb*DSIx(:,par.pindx.lbeta)] + ...
              [-d0(Gx(:,par.pindx.lbeta))*SILx(:,par.pindx.lbsi).*Si2C; ...
               d0(Gx(:,par.pindx.lbeta))*SILx(:,par.pindx.lbsi).*Si2C];
        Sixx(:,kk) = mfactor(FD, tmp);
        kk = kk + 1;
    end
    % beta at
    if (par.opt_beta == on & par.opt_at == on)
        tmp = at*[d0(dkdat)*DSIx(:,par.pindx.lbeta); ...
                  -(d0(dkdat)+dPFDdat)*DSIx(:,par.pindx.lbeta)] + ...
              [-d0(Gx(:,par.pindx.lbeta))*SILx(:,par.pindx.lat).*Si2C; ...
               d0(Gx(:,par.pindx.lbeta))*SILx(:,par.pindx.lat).*Si2C];

        Sixx(:,kk) = mfactor(FD, tmp);
        kk = kk + 1;
    end
    % beta bt
    if (par.opt_beta == on & par.opt_bt == on)
        tmp = bt*[d0(dkdbt)*DSIx(:,par.pindx.lbeta); ...
                  -(d0(dkdbt)+dPFDdbt)*DSIx(:,par.pindx.lbeta)] + ...
              [-d0(Gx(:,par.pindx.lbeta))*SILx(:,par.pindx.lbt).*Si2C; ...
               d0(Gx(:,par.pindx.lbeta))*SILx(:,par.pindx.lbt).*Si2C];

        Sixx(:,kk) = mfactor(FD, tmp);
        kk = kk + 1;
    end
    % beta aa
    if (par.opt_beta == on & par.opt_aa == on)
        tmp = [-d0(Gxx(:,kk))*SIL; ...
               d0(Gxx(:,kk))*SIL] + ...
              [-d0(G*SILx(:,par.pindx.lbeta))*dSi2Cdaa; ...
               d0(G*SILx(:,par.pindx.lbeta))*dSi2Cdaa] + ...
              [-d0(Gx(:,par.pindx.lbeta))*SILx(:,par.pindx.aa).*Si2C; ...
               d0(Gx(:,par.pindx.lbeta))*SILx(:,par.pindx.aa).*Si2C];

        Sixx(:,kk) = mfactor(FD, tmp);
        kk = kk + 1;
    end
    % beta bb
    if (par.opt_beta == on & par.opt_bb == on)
        tmp = [-d0(Gxx(:,kk))*SIL; ...
               d0(Gxx(:,kk))*SIL] + ...
              bb*[-d0(G*SILx(:,par.pindx.lbeta))*dSi2Cdbb; ...
                  d0(G*SILx(:,par.pindx.lbeta))*dSi2Cdbb] + ...
              [-d0(Gx(:,par.pindx.lbeta))*SILx(:,par.pindx.lbb).*Si2C; ...
               d0(Gx(:,par.pindx.lbeta))*SILx(:,par.pindx.lbb).*Si2C];

        Sixx(:,kk) = mfactor(FD, tmp);
        kk = kk + 1;
    end
    % bsi bsi
    if (par.opt_bsi == on)
        [~,~,hout] = buildPFD(M3d,grd,par,'bSi');
        d2PFDdb2 = hout.PFD_b_b;
        tmp = bsi*[0*DSI; ...
                   -dPFDdb*DSI] + ...
              bsi*bsi*[0*DSI; ...
                       -d2PFDdb2*DSI] + ...
              2*bsi*[0*DSI; ...
                     -dPFDdb*DSIx(:,par.pindx.lbsi)];

        Sixx(:,kk) = mfactor(FD, tmp);
        kk = kk + 1;
    end
    % bsi at
    if (par.opt_bsi & par.opt_at)
        [~,~,hout] = buildPFD(M3d,grd,par,'bSi');
        d2PFDdatdb = hout.PFD_at_b;
        tmp = bsi*at*[0*DSI; ...
                      -d2PFDdatdb*DSI] + ...
              bsi*[0*DSI; ...
                   -dPFDdb*DSIx(:,par.pindx.lat)] + ...
              at*[d0(dkdat)*DSIx(:,par.pindx.lbsi); ...
                  -(d0(dkdat)+dPFDdat)*DSIx(:,par.pindx.lbsi)];
        
        Sixx(:,kk) = mfactor(FD, tmp);
        kk = kk + 1;
    end
    % bsi bt
    if (par.opt_bsi & par.opt_bt)
        vout = buildPFD(par, par);
        d2PFDdbtdb = vout.d2PFDdbtdb;
        tmp = bsi*bt*[0*DSI; ...
                      -d2PFDdbtdb*DSI] + ...
              bsi*[0*DSI; ...
                   -dPFDdb*DSIx(:,par.pindx.lbt)] + ...
              bt*[d0(dkdbt)*DSIx(:,par.pindx.lbsi); ...
                  -(d0(dkdbt)+dPFDdbt)*DSIx(:,par.pindx.lbsi)];
        
        Sixx(:,kk) = mfactor(FD, tmp);
        kk = kk + 1;
    end
    % bsi aa
    if (par.opt_bsi & par.opt_aa)
        tmp = bsi*[0*DSI; ...
                   -dPFDdb*DSIx(:,par.pindx.aa)] + ...
              [-d0(G*SILx(:,par.pindx.lbsi))*dSi2Cdaa; ...
               d0(G*SILx(:,par.pindx.lbsi))*dSi2Cdaa] + ...
              [-d0(Gx(:,par.pindx.lbsi))*SILx(:,par.pindx.aa).*Si2C; ...
               d0(Gx(:,par.pindx.lbsi))*SILx(:,par.pindx.aa).*Si2C];
        
        Sixx(:,kk) = mfactor(FD, tmp);
        kk = kk + 1;
    end
    % bsi bb
    if (par.opt_bsi & par.opt_bb)
        tmp = bsi*[0*DSI; ...
                   -dPFDdb*DSIx(:,par.pindx.lbb)] + ...
              bb*[-d0(G*SILx(:,par.pindx.lbsi))*dSi2Cdbb; ...
                  d0(G*SILx(:,par.pindx.lbsi))*dSi2Cdbb] + ...
              [-d0(Gx(:,par.pindx.lbsi))*SILx(:,par.pindx.lbb).*Si2C; ...
               d0(Gx(:,par.pindx.lbsi))*SILx(:,par.pindx.lbb).*Si2C];

        Sixx(:,kk) = mfactor(FD, tmp);
        kk = kk + 1;
    end
    % at at
    if (par.opt_at == on)
        [~,~,hout] = buildPFD(M3d,grd,par,'bSi');
        d2PFDdat2 = hout.PFD_at_at;
        d2kdat2 = dkdat;
        tmp = at*[d0(dkdat)*DSI; ...
                  -(d0(dkdat)+dPFDdat)*DSI] + ...
              at*at*[0*DSI; ...
                     -d2PFDdat2*DSI] + ...
              2*at*[d0(dkdat)*DSIx(:,par.pindx.lat); ...
                  -(d0(dkdat)+dPFDdat)*DSIx(:,par.pindx.lat)];

        Sixx(:,kk) = mfactor(FD, tmp);
        kk = kk + 1;
    end
    % at bt
    if (par.opt_at == on & par.opt_bt == on)
        [~,~,hout] = buildPFD(M3d,grd,par,'bSi');
        d2PFDdatdbt = hout.PFD_at_bt;
        d2kdatdbt = -exp(-bt./T).*(1./T);
        tmp = at*bt*[d0(d2kdatdbt); ...
                     -(d0(d2kdatdbt)+d2PFDdatdbt)]*DSI + ...
              at*[d0(dkdat)*DSIx(:,par.pindx.lbt); ...
                  -(d0(dkdat)+dPFDdat)*DSIx(:,par.pindx.lbt)]+ ...
              bt*[d0(dkdbt)*DSIx(:,par.pindx.lat); ...
                  -(d0(dkdbt)+dPFDdbt)*DSIx(:,par.pindx.lat)];
        
        Sixx(:,kk) = mfactor(FD, tmp);
        kk = kk + 1;
    end
    % at aa
    if (par.opt_at == on & par.opt_aa == on)
        tmp = at*[d0(dkdat)*DSIx(:,par.pindx.aa); ...
                  -(d0(dkdat)+dPFDdat)*DSIx(:,par.pindx.aa)] + ...
              [-d0(G*SILx(:,par.pindx.lat))*dSi2Cdaa; ... 
               d0(G*SILx(:,par.pindx.lat))*dSi2Cdaa] + ...
              [-d0(Gx(:,par.pindx.lat))*SILx(:,par.pindx.aa).*Si2C; ...
               d0(Gx(:,par.pindx.lat))*SILx(:,par.pindx.aa).*Si2C];
        
        Sixx(:,kk) = mfactor(FD, tmp);
        kk = kk + 1;
    end
    % at bb
    if (par.opt_at == on & par.opt_bb == on)
        tmp = at*[d0(dkdat)*DSIx(:,par.pindx.lbb); ...
                  -(d0(dkdat)+dPFDdat)*DSIx(:,par.pindx.lbb)] + ...
              bb*[-d0(G*SILx(:,par.pindx.lat))*dSi2Cdbb; ... 
                  d0(G*SILx(:,par.pindx.lat))*dSi2Cdbb] + ...
              [-d0(Gx(:,par.pindx.lat))*SILx(:,par.pindx.lbb).*Si2C; ...
               d0(Gx(:,par.pindx.lat))*SILx(:,par.pindx.lbb).*Si2C];

        Sixx(:,kk) = mfactor(FD, tmp);
        kk = kk + 1;
    end
    % bt bt
    if (par.opt_bt == on)
        [~,~,hout] = buildPFD(M3d,grd,par,'bSi');
        d2kdbtdbt = kappa_si./T.^2;
        d2PFDdbt2 = hout.PFD_bt_bt;
        tmp = bt*[d0(dkdbt)*DSI; ...
                  -(d0(dkdbt)+dPFDdbt)*DSI] + ...
              bt*bt*[d0(d2kdbtdbt)*DSI; ...
                     -(d0(d2kdbtdbt)+d2PFDdbt2)*DSI] + ...
              2*bt*[d0(dkdbt)*DSIx(:,par.pindx.lbt); ...
                    -(d0(dkdbt)+dPFDdbt)*DSIx(:,par.pindx.lbt)];
        
        Sixx(:,kk) = mfactor(FD, tmp);
        kk = kk + 1;
    end
    % bt aa
    if (par.opt_bt == on & par.opt_aa == on)
        tmp = bt*[d0(dkdbt)*DSIx(:,par.pindx.aa); ...
                  -(d0(dkdbt) + dPFDdbt)*DSIx(:,par.pindx.aa)] + ...
              [-d0(G*SILx(:,par.pindx.lbt))*dSi2Cdaa; ... 
               d0(G*SILx(:,par.pindx.lbt))*dSi2Cdaa] + ...
              [-d0(Gx(:,par.pindx.lbt))*SILx(:,par.pindx.aa).*Si2C; ...
               d0(Gx(:,par.pindx.lbt))*SILx(:,par.pindx.aa).*Si2C];
        
        Sixx(:,kk) = mfactor(FD, tmp);
        kk = kk + 1;
    end
    % bt bb
    if (par.opt_bt == on & par.opt_bb == on)
        tmp = bt*[d0(dkdbt)*DSIx(:,par.pindx.lbb); ...
                  -(d0(dkdbt)+dPFDdbt)*DSIx(:,par.pindx.lbb)] + ...
              bb*[-d0(G*SILx(:,par.pindx.lbt))*dSi2Cdbb; ... 
                  d0(G*SILx(:,par.pindx.lbt))*dSi2Cdbb] + ...
              [-d0(Gx(:,par.pindx.lbt))*SILx(:,par.pindx.lbb).*Si2C; ...
               d0(Gx(:,par.pindx.lbt))*SILx(:,par.pindx.lbb).*Si2C];

        Sixx(:,kk) = mfactor(FD, tmp);
        kk = kk + 1;
    end
    % aa aa
    if (par.opt_aa == on & par.opt_aa == on)
        tmp = [-d0(Gxx(:,kk))*SIL; ...
               d0(Gxx(:,kk))*SIL] + ...
              2*[-d0(G*SILx(:,par.pindx.aa))*dSi2Cdaa; ...
                 d0(G*SILx(:,par.pindx.aa))*dSi2Cdaa];
        
        Sixx(:,kk) = mfactor(FD, tmp);
        kk = kk + 1;
    end
    % aa bb
    if (par.opt_aa == on & par.opt_bb == on)
        tmp = [-d0(Gxx(:,kk))*SIL; ...
               d0(Gxx(:,kk))*SIL] + ...
              bb*[-d0(G*SILx(:,par.pindx.aa))*dSi2Cdbb; ... 
                  d0(G*SILx(:,par.pindx.aa))*dSi2Cdbb] + ...
              [-d0(G*SILx(:,par.pindx.lbb))*dSi2Cdaa; ... 
               d0(G*SILx(:,par.pindx.lbb))*dSi2Cdaa];
              
        Sixx(:,kk) = mfactor(FD, tmp);
        kk = kk + 1;
    end
    % bb bb
    if (par.opt_bb == on & par.opt_bb == on)
        tmp = [-d0(Gxx(:,kk))*SIL; ...
               d0(Gxx(:,kk))*SIL] + ...
              2*bb*[-d0(G*SILx(:,par.pindx.lbb))*dSi2Cdbb; ...
                    d0(G*SILx(:,par.pindx.lbb))*dSi2Cdbb];
        
        Sixx(:,kk) = mfactor(FD, tmp);
        kk = kk + 1;
    end
end

end

