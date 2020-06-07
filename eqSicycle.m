function [parm,Si,Six,Sixx] = eqSicycle(par, parm, x)
on = true; off = false;
M3d =  parm.M3d;
iwet = parm.iwet;
nwet = parm.nwet;
I    = speye(nwet); % make an identity matrix;
parm.nx = length(x);
nsx = 0; % count number of Si parameters
% unpack the parameters to be optimized
if (par.opt_alpha == on)
    lalpha = x(par.pindx.lalpha);
    parm.alpha  = exp(lalpha);
else
    parm.alpha  = par.alpha;
end
% beta
if (par.opt_beta == on)
    lbeta = x(par.pindx.lbeta);
    parm.beta  = exp(lbeta);
else
    parm.beta  = par.beta;
end
% bsi
if (par.opt_bsi == on)
    nsx = nsx + 1;
    lbsi = x(par.pindx.lbsi);
    parm.bsi = exp(lbsi);
else
    parm.bsi  = par.bsi;
end
% at
if (par.opt_at == on)
    nsx = nsx + 1;
    lat = x(par.pindx.lat);
    parm.at = exp(lat);
else
    parm.at  = par.at;
end
% bt
if (par.opt_bt == on)
    nsx = nsx + 1;
    lbt = x(par.pindx.lbt);
    parm.bt  = exp(lbt);
else
    parm.bt  = par.bt;
end
% aa
if (par.opt_aa == on)
    nsx = nsx + 1;
    parm.aa = x(par.pindx.aa);
else
    parm.aa  = par.aa;
end
% bb
if (par.opt_bb == on)
    nsx = nsx + 1;
    lbb = x(par.pindx.lbb);
    parm.bb  = exp(lbb);
else
    parm.bb  = par.bb;
end
% bb
if (par.opt_kappa_gs == on)
    nsx = nsx + 1;
    lkappa_gs = x(par.pindx.lkappa_gs);
    parm.kappa_gs  = exp(lkappa_gs);
else
    parm.kappa_gs  = par.kappa_gs;
end
parm.nsx = nsx;
%
% ++++++++++++++++++++++++++++++++++++++++++++++++++
at  = parm.at;
bt  = parm.bt;
aa  = parm.aa;
bb  = parm.bb;
bsi = parm.bsi;
kappa_gs = parm.kappa_gs;
T   = parm.sst(iwet) + 273.15;
kappa_si = at * exp(-bt./T); 
% ++++++++++++++++++++++++++++++++++++++++++++++++++
%
% fixed parameters
SI4 = parm.SIL(iwet);
SILbar  = parm.SILbar*M3d(iwet);
%
%%%%%%%%%%%%%%%%%%% create Si2P %%%%%%%%%%%%%%%%%%%%%
smsk = M3d;
smsk(:,:,3:end) = 0;
isrf = find(smsk(iwet));
dVs = parm.dVt(iwet(isrf));
surface_mean = @(x) sum(x(isrf).*dVs)/sum(dVs);

% compute the mean of the regressor variable
Z = SI4;
mu = surface_mean(Z);
Delta = sqrt(surface_mean((Z-mu).^2));

% standardize the regressor variables
ZR = 0.5+0.5*tanh((Z-mu)/Delta);
% ZR = (Z-min(Z(isrf)))./(max(Z(isrf))-min(Z(isrf)));
%
Si2C = (aa*ZR + bb);      parm.Si2C = Si2C;
dSi2Cdaa = ZR;            parm.dSi2Cdaa = dSi2Cdaa;
dSi2Cdbb = ones(nwet,1);  parm.dSi2Cdbb = dSi2Cdbb;
%%%%%%%%%%%%%%%%%%%%%%% end %%%%%%%%%%%%%%%%%%%%%%%%
%
vout  = buildPFD(par,parm);
TRdiv = parm.TRdiv;
PFdiv = vout.PFdiv;
G     = uptake(par, parm);
% build Jacobian matrix
% +++++++++++++++++++++++++++++++++++++++++
% F = ([TRdiv*SIL + G*Si2C - kappa_si*DSI + kappa_gs*(SIL-SILbar);...
% PFdiv*DSI - G*Si2C + kappa_si*DSI]);
Jac = [[TRdiv+kappa_gs*I, -d0(kappa_si)];...
       [0*I,  PFdiv+d0(kappa_si)]];

RHS = [[SILbar*kappa_gs - G.*Si2C]; ...
       [sparse(nwet,1) + G.*Si2C]];

FD = mfactor(Jac);
Si = mfactor(FD,RHS);
SIL = Si(1:nwet);
DSI = Si(nwet+1:end);
save tmpSi SIL DSI
%% +++++++++++++++++++++++++++++++++++++++++++
if nargout > 1
    [~,Gx] = uptake(par, parm);
    % sigma
    if (par.opt_sigma == on)
        tmp =  -[Si2C.*Gx(:,par.pindx.lsigma) ;...
                -Si2C.*Gx(:,par.pindx.lsigma)]   ;
        Six(:,par.pindx.lsigma) = mfactor(FD, tmp);
    end
    
    % interpp
    if (par.opt_interpp == on)
        tmp = -[Si2C.*Gx(:,par.pindx.linterpp) ;...
                -Si2C.*Gx(:,par.pindx.linterpp)]   ;
        Six(:,par.pindx.linterpp) = mfactor(FD, tmp);
    end
    
    % slopep
    if (par.opt_slopep == on)
        tmp = -[Si2C.*Gx(:,par.pindx.slopep) ;...
                -Si2C.*Gx(:,par.pindx.slopep)]   ;
        Six(:,par.pindx.slopep) = mfactor(FD, tmp);
    end
    
    % kappa_dp
    if (par.opt_kappa_dp == on)
        tmp = -[Si2C.*Gx(:,par.pindx.lkappa_dp) ;...
                -Si2C.*Gx(:,par.pindx.lkappa_dp)]   ;
        Six(:,par.pindx.lkappa_dp) = mfactor(FD, tmp);
    end

    % alpha
    if (par.opt_alpha == on)
        tmp = -[Si2C.*Gx(:,par.pindx.lalpha); ...
                -Si2C.*Gx(:,par.pindx.lalpha)];
        
        Six(:,par.pindx.lalpha) = mfactor(FD, tmp);
    end

    % beta
    if (par.opt_beta == on)
        tmp = -[Si2C.*Gx(:,par.pindx.lbeta); ...
                -Si2C.*Gx(:,par.pindx.lbeta)];

        Six(:,par.pindx.lbeta) = mfactor(FD, tmp);
    end
    % bsi
    if (par.opt_bsi == on)
        vout = buildPFD(par, parm);
        dPFDdb = vout.dPFDdb;
        tmp = -bsi*[DSI*0; ...
                    dPFDdb*DSI];
        
        Six(:,par.pindx.lbsi) = mfactor(FD, tmp);
    end
    % at
    if (par.opt_at == on)
        vout = buildPFD(par, parm);
        dPFDdat = vout.dPFDdat;
        dkdat = exp(-bt./T);
        tmp = -at*[-dkdat.*DSI; ...
                   dkdat.*DSI + dPFDdat*DSI];
        
        Six(:,par.pindx.lat) = mfactor(FD, tmp);
    end
    % bt
    if (par.opt_bt == on)
        vout = buildPFD(par, parm);
        dPFDdbt = vout.dPFDdbt;
        dkdbt = -(at*exp(-bt./T))./T;
        tmp = -bt*[-dkdbt.*DSI; ...
                   dkdbt.*DSI + dPFDdbt*DSI];

        Six(:,par.pindx.lbt) = mfactor(FD, tmp);
    end
    % aa
    if (par.opt_aa == on)
        tmp = -[d0(G)*dSi2Cdaa; ... %.*SIL; ...
                -d0(G)*dSi2Cdaa]; %.*SIL];

        Six(:,par.pindx.aa) = mfactor(FD, tmp);
    end
    % bb
    if (par.opt_bb == on)
        tmp = -bb*[d0(G)*dSi2Cdbb; ... %.*SIL; ...
                   -d0(G)*dSi2Cdbb]; %.*SIL];

        Six(:,par.pindx.lbb) = mfactor(FD, tmp);
    end
    % kappa_gs
    if (par.opt_kappa_gs == on)
        tmp = kappa_gs*[SILbar-SIL; ...
                        sparse(nwet,1)];

        Six(:,par.pindx.lkappa_gs) = mfactor(FD, tmp);
    end    
end
SILx = Six(1:nwet,:); DSIx = Six(nwet+1:end,:);
%% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
if (nargout > 2)
    nx = parm.npx + parm.nsx;
    ncs = nchoosek(nx,2)+nx;
    Sixx = sparse(2*nwet,ncs);
    % Compute the hessian of the solution wrt the parameters
    DIPxx = parm.Pxx(1:nwet,:);
    [~,~,Gxx] = uptake(par, parm);
    % sigma sigma
    kk = 1;
    if (par.opt_sigma == on)
        tmp = [-d0(Gxx(:,kk))*Si2C; ...
               d0(Gxx(:,kk))*Si2C];
        Sixx(:,kk) = mfactor(FD, tmp);
        kk = kk + 1;
    end
    % sigma kappa_dp
    if (par.opt_sigma == on & par.opt_kappa_dp == on)
        tmp = [-d0(Gxx(:,kk))*Si2C; ...
               d0(Gxx(:,kk))*Si2C];
        Sixx(:,kk) = mfactor(FD, tmp);
        kk = kk + 1;
    end
    % sigma slopep
    if (par.opt_sigma == on & par.opt_slopep == on)
        tmp = [-d0(Gxx(:,kk))*Si2C; ...
               d0(Gxx(:,kk))*Si2C];
        Sixx(:,kk) = mfactor(FD, tmp);
        kk = kk + 1;
    end
    % sigma interpp
    if (par.opt_sigma == on & par.opt_interpp == on)
        tmp = [-d0(Gxx(:,kk))*Si2C; ...
               d0(Gxx(:,kk))*Si2C];
        Sixx(:,kk) = mfactor(FD, tmp);
        kk = kk + 1;
    end
    % sigma alpha
    if (par.opt_sigma == on & par.opt_alpha == on)
        tmp = [-d0(Gxx(:,kk))*Si2C; ...
               d0(Gxx(:,kk))*Si2C];
        Sixx(:,kk) = mfactor(FD, tmp);
        kk = kk + 1;
    end
    % sigma beta
    if (par.opt_sigma == on & par.opt_beta == on)
        tmp = [-d0(Gxx(:,kk))*Si2C; ...
               d0(Gxx(:,kk))*Si2C];
        Sixx(:,kk) = mfactor(FD, tmp);
        kk = kk + 1;
    end
    % kappa_dp kappa_dp
    if (par.opt_kappa_dp == on)
        tmp = [-d0(Gxx(:,kk))*Si2C; ...
               d0(Gxx(:,kk))*Si2C];
        Sixx(:,kk) = mfactor(FD, tmp);
        kk = kk + 1;
    end
    % kappa_dp slopep
    if (par.opt_kappa_dp == on & par.opt_slopep == on)
        tmp = [-d0(Gxx(:,kk))*Si2C; ...
               d0(Gxx(:,kk))*Si2C];
        Sixx(:,kk) = mfactor(FD, tmp);
        kk = kk + 1;
    end
    % kappa_dp interpp
    if (par.opt_kappa_dp == on & par.opt_interpp == on)
        tmp = [-d0(Gxx(:,kk))*Si2C; ...
               d0(Gxx(:,kk))*Si2C];
        Sixx(:,kk) = mfactor(FD, tmp);
        kk = kk + 1;
    end
    % kappa_dp alpha
    if (par.opt_kappa_dp == on & par.opt_alpha == on)
        tmp = [-d0(Gxx(:,kk))*Si2C; ...
               d0(Gxx(:,kk))*Si2C];
        Sixx(:,kk) = mfactor(FD, tmp);
        kk = kk + 1;
    end
    % kappa_dp beta
    if (par.opt_kappa_dp == on & par.opt_beta == on)
        tmp = [-d0(Gxx(:,kk))*Si2C; ...
               d0(Gxx(:,kk))*Si2C];
        Sixx(:,kk) = mfactor(FD, tmp);
        kk = kk + 1;
    end
    % slopep slopep
    if (par.opt_slopep == on)
        tmp = [-d0(Gxx(:,kk))*Si2C; ...
               d0(Gxx(:,kk))*Si2C];
        Sixx(:,kk) = mfactor(FD, tmp);
        kk = kk + 1;
    end
    % slopep interpp
    if (par.opt_slopep == on & par.opt_interpp == on)
        tmp = [-d0(Gxx(:,kk))*Si2C; ...
               d0(Gxx(:,kk))*Si2C];
        Sixx(:,kk) = mfactor(FD, tmp);
        kk = kk + 1;
    end
    % slopep alpha
    if (par.opt_slopep == on & par.opt_alpha == on)
        tmp = [-d0(Gxx(:,kk))*Si2C; ...
               d0(Gxx(:,kk))*Si2C];
        Sixx(:,kk) = mfactor(FD, tmp);
        kk = kk + 1;
    end
    % slopep beta
    if (par.opt_slopep == on & par.opt_beta == on)
        tmp = [-d0(Gxx(:,kk))*Si2C; ...
               d0(Gxx(:,kk))*Si2C];
        Sixx(:,kk) = mfactor(FD, tmp);
        kk = kk + 1;
    end
    % interpp interpp
    if (par.opt_interpp == on)
        tmp = [-d0(Gxx(:,kk))*Si2C; ...
               d0(Gxx(:,kk))*Si2C];
        Sixx(:,kk) = mfactor(FD, tmp);
        kk = kk + 1;
    end
    % interpp alpha
    if (par.opt_interpp == on & par.opt_alpha == on)
        tmp = [-d0(Gxx(:,kk))*Si2C; ...
               d0(Gxx(:,kk))*Si2C];
        Sixx(:,kk) = mfactor(FD, tmp);
        kk = kk + 1;
    end
    % interpp beta
    if (par.opt_interpp == on & par.opt_beta == on)
        tmp = [-d0(Gxx(:,kk))*Si2C; ...
               d0(Gxx(:,kk))*Si2C];
        Sixx(:,kk) = mfactor(FD, tmp);
        kk = kk + 1;
    end
    % alpha alpha
    if (par.opt_alpha == on)
        tmp = [-d0(Si2C)*Gxx(:,kk); ...
               d0(Si2C)*Gxx(:,kk)]; %drhsdalphadalpha      
        Sixx(:,kk) = mfactor(FD, tmp);
        kk = kk + 1;
    end
    % alpha beta
    if (par.opt_alpha == on & par.opt_beta == on)
        tmp = [-d0(Si2C)*Gxx(:,kk); ...
               d0(Si2C)*Gxx(:,kk)]; %drhsdalphadbeta
        Sixx(:,kk) = mfactor(FD, tmp);
        kk = kk + 1;
    end
    % beta beta
    if (par.opt_beta == on)
        tmp = [-d0(Si2C)*Gxx(:,kk); ...
               d0(Si2C)*Gxx(:,kk)]; %drhsdalphadbeta
        Sixx(:,kk) = mfactor(FD, tmp);
        kk = kk + 1;
    end
    % --------------------------------------------------------
    % Simodel parameters
    % sigma bsi
    if (par.opt_sigma == on & par.opt_bsi == on)
        tmp = bsi*[0*DSI; ...
                  -dPFDdb*DSIx(:,par.pindx.lsigma)];
        Sixx(:,kk) = mfactor(FD, tmp);
        kk = kk + 1;
    end
    % sigma at
    if (par.opt_sigma == on & par.opt_at == on)
        tmp = at*[d0(dkdat)*DSIx(:,par.pindx.lsigma); ...
                  -(d0(dkdat) + dPFDdat)*DSIx(:,par.pindx.lsigma)];
        Sixx(:,kk) = mfactor(FD, tmp);
        kk = kk + 1;
    end
    % sigma bt
    if (par.opt_sigma == on & par.opt_bt == on)
        tmp = bt*[d0(dkdbt)*DSIx(:,par.pindx.lsigma); ...
                  -(d0(dkdbt)+dPFDdbt)*DSIx(:,par.pindx.lsigma)];
        Sixx(:,kk) = mfactor(FD, tmp);
        kk = kk + 1;
    end
    % sigma aa
    if (par.opt_sigma == on & par.opt_aa == on)
        tmp = [-Gxx(:,kk); ...
               Gxx(:,kk)];
        Sixx(:,kk) = mfactor(FD, tmp);
        kk = kk + 1;
    end
    % sigma bb
    if (par.opt_sigma == on & par.opt_bb == on)
        tmp = [-Gxx(:,kk); ...
               Gxx(:,kk)];
        Sixx(:,kk) = mfactor(FD, tmp);
        kk = kk + 1;
    end
    % kappa_dp bsi
    if (par.opt_kappa_dp == on & par.opt_bsi == on)
        tmp = bsi*[0*DSI; ...
                  -(dPFDdb)*DSIx(:,par.pindx.lkappa_dp)];
        Sixx(:,kk) = mfactor(FD, tmp);
        kk = kk + 1;
    end
    % kappa_dp at
    if (par.opt_kappa_dp == on & par.opt_at == on)
        tmp = at*[d0(dkdat)*DSIx(:,par.pindx.lkappa_dp); ...
                  -(d0(dkdat)+dPFDdat)*DSIx(:,par.pindx.lkappa_dp)];
        Sixx(:,kk) = mfactor(FD, tmp);
        kk = kk + 1;
    end
    % kappa_dp bt
    if (par.opt_kappa_dp == on & par.opt_bt == on)
        tmp = bt*[d0(dkdbt)*DSIx(:,par.pindx.lkappa_dp); ...
                  -(d0(dkdbt)+dPFDdbt)*DSIx(:,par.pindx.lkappa_dp)];
        Sixx(:,kk) = mfactor(FD, tmp);
        kk = kk + 1;
    end
    % kappa_dp aa
    if (par.opt_kappa_dp == on & par.opt_aa == on)
        tmp =  [-Gxx(:,kk); ...
                Gxx(:,kk)];
        Sixx(:,kk) = mfactor(FD, tmp);
        kk = kk + 1;
    end
    % kappa_dp bb
    if (par.opt_kappa_dp == on & par.opt_bb == on)
        tmp = [-Gxx(:,kk); ...
               Gxx(:,kk)];
        Sixx(:,kk) = mfactor(FD, tmp);
        kk = kk + 1;
    end
    % slope bsi
    if (par.opt_slopep == on & par.opt_bsi == on)
        tmp = bsi*[0*DSI; ...
                   -dPFDdb*DSIx(:,par.pindx.slopep)];
        Sixx(:,kk) = mfactor(FD, tmp);
        kk = kk + 1;
    end
    % slope at
    if (par.opt_slopep == on & par.opt_at == on)
        tmp = at*[d0(dkdat)*DSIx(:,par.pindx.slopep); ...
                  -(d0(dkdat)+dPFDdat)*DSIx(:,par.pindx.slopep)];
        Sixx(:,kk) = mfactor(FD, tmp);
        kk = kk + 1;
    end
    % slope bt
    if (par.opt_slopep == on & par.opt_bt == on)
        tmp = bt*[d0(dkdbt)*DSIx(:,par.pindx.slopep); ...
                  -(d0(dkdbt)+dPFDdbt)*DSIx(:,par.pindx.slopep)];
        Sixx(:,kk) = mfactor(FD, tmp);
        kk = kk + 1;
    end
    % slopep aa
    if (par.opt_slopep == on & par.opt_aa == on)
        tmp = [-Gxx(:,kk); ...
               Gxx(:,kk)];
        Sixx(:,kk) = mfactor(FD, tmp);
        kk = kk + 1;
    end
    % slopep bb
    if (par.opt_slopep == on & par.opt_bb == on)
        tmp = [-Gxx(:,kk); ...
               Gxx(:,kk)];
        Sixx(:,kk) = mfactor(FD, tmp);
        kk = kk + 1;
    end
    % interpp bsi
    if (par.opt_interpp == on & par.opt_bsi == on)
        tmp = bsi*[0*DSI; ...
                   -dPFDdb*DSIx(:, par.pindx.linterpp)];
        Sixx(:,kk) = mfactor(FD, tmp);
        kk = kk + 1;
    end
    % interpp at
    if (par.opt_interpp == on & par.opt_at == on)
        tmp = at*[d0(dkdat)*DSIx(:, par.pindx.linterpp); ...
                  -(d0(dkdat)+dPFDdat)*DSIx(:, par.pindx.linterpp)];
        Sixx(:,kk) = mfactor(FD, tmp);
        kk = kk + 1;
    end
    % interpp bt
    if (par.opt_interpp == on & par.opt_bt == on)
        tmp = bt*[d0(dkdbt)*DSIx(:, par.pindx.linterpp); ...
                  -(d0(dkdbt)+dPFDdbt)*DSIx(:,par.pindx.linterpp)];
        Sixx(:,kk) = mfactor(FD, tmp);
        kk = kk + 1;
    end
    % interpp aa
    if (par.opt_interpp == on & par.opt_aa == on)
        tmp = [-Gxx(:,kk); ...
               Gxx(:,kk)];
        Sixx(:,kk) = mfactor(FD, tmp);
        kk = kk + 1;
    end
    % interpp bb
    if (par.opt_interpp == on & par.opt_bb == on)
        tmp = [-Gxx(:,kk); ...
               Gxx(:,kk)];
        Sixx(:,kk) = mfactor(FD, tmp);
        kk = kk + 1;
    end
    % alpha bsi
    if (par.opt_alpha == on & par.opt_bsi == on)
        tmp = bsi*[0*DSI; ...
                   -dPFDdb*DSIx(:,par.pindx.lalpha)];
        Sixx(:,kk) = mfactor(FD, tmp);
        kk = kk + 1;
    end
    % alpha at
    if (par.opt_alpha == on & par.opt_at == on)
        tmp = at*[d0(dkdat)*DSIx(:,par.pindx.lalpha); ...
                  -(d0(dkdat)+dPFDdat)*DSIx(:,par.pindx.lalpha)];
        Sixx(:,kk) = mfactor(FD, tmp);
        kk = kk + 1;
    end
    % alpha bt
    if (par.opt_alpha == on & par.opt_bt == on)
        tmp = bt*[d0(dkdbt)*DSIx(:,par.pindx.lalpha); ...
                  -(d0(dkdbt) + dPFDdbt)*DSIx(:,par.pindx.lalpha)];
        Sixx(:,kk) = mfactor(FD, tmp);
        kk = kk + 1;
    end
    % alpha aa
    if (par.opt_alpha == on & par.opt_aa == on)
        tmp = [-Gxx(:,kk); ...
               Gxx(:,kk)];
        Sixx(:,kk) = mfactor(FD, tmp);
        kk = kk + 1;
    end
    % alpha bb
    if (par.opt_alpha == on & par.opt_bb == on)
        tmp = [-Gxx(:,kk); ...
               Gxx(:,kk)];
        Sixx(:,kk) = mfactor(FD, tmp);
        kk = kk + 1;
    end
    % beta bsi
    if (par.opt_beta == on & par.opt_bsi == on)
        tmp = bsi*[0*DSI; ...
                   -dPFDdb*DSIx(:,par.pindx.lbeta)];
        Sixx(:,kk) = mfactor(FD, tmp);
        kk = kk + 1;
    end
    % beta at
    if (par.opt_beta == on & par.opt_at == on)
        tmp = at*[d0(dkdat)*DSIx(:,par.pindx.lbeta); ...
                  -(d0(dkdat)+dPFDdat)*DSIx(:,par.pindx.lbeta)];
        Sixx(:,kk) = mfactor(FD, tmp);
        kk = kk + 1;
    end
    % beta bt
    if (par.opt_beta == on & par.opt_bt == on)
        tmp = bt*[d0(dkdbt)*DSIx(:,par.pindx.lbeta); ...
                  -(d0(dkdbt)+dPFDdbt)*DSIx(:,par.pindx.lbeta)];
        Sixx(:,kk) = mfactor(FD, tmp);
        kk = kk + 1;
    end
    % beta aa
    if (par.opt_beta == on & par.opt_aa == on)
        tmp = [-Gxx(:,kk); ...
               Gxx(:,kk)];
        Sixx(:,kk) = mfactor(FD, tmp);
        kk = kk + 1;
    end
    % beta bb
    if (par.opt_beta == on & par.opt_bb == on)
        tmp = [-Gxx(:,kk); ...
               Gxx(:,kk)];
        Sixx(:,kk) = mfactor(FD, tmp);
        kk = kk + 1;
    end
    % bsi bsi
    if (par.opt_bsi == on)
        vout = buildPFD(par, parm);
        d2PFDdb2 = vout.d2PFDdb2;
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
        vout = buildPFD(par, parm);
        d2PFDdatdb = vout.d2PFDdatdb;
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
        vout = buildPFD(par, parm);
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
              [-d0(Gx(:,par.pindx.lbsi))*dSi2Cdaa; ... %.*SIL; ...
               d0(Gx(:,par.pindx.lbsi))*dSi2Cdaa]; %.*SIL];
        
        Sixx(:,kk) = mfactor(FD, tmp);
        kk = kk + 1;
    end
    % bsi bb
    if (par.opt_bsi & par.opt_bb)
        tmp = bsi*[0*DSI; ...
                   -dPFDdb*DSIx(:,par.pindx.lbb)] + ...
              [-d0(Gx(:,par.pindx.lbsi))*dSi2Cdbb; ... %.*SIL; ...
               d0(Gx(:,par.pindx.lbsi))*dSi2Cdbb]; %.*SIL];

        Sixx(:,kk) = mfactor(FD, tmp);
        kk = kk + 1;
    end
    % at at
    if (par.opt_at == on)
        vout = buildPFD(par, parm);
        d2PFDdat2 = vout.d2PFDdat2;
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
        vout = buildPFD(par, parm);
        d2PFDdatdbt = vout.d2PFDdatdbt;
        d2kdatdbt = -exp(-bt./T).*(1./T);
        tmp = at*bt*[d0(d2kdatdbt); ...
                     -(d0(d2kdatdbt)+d2PFDdatdbt)]*DSI + ...
              at*[d0(dkdat)*DSIx(:,par.pindx.lbt); ...
                  -(d0(dkdat)+dPFDdat)*DSIx(:,par.pindx.lbt)] + ...
              bt*[d0(dkdbt)*DSIx(:,par.pindx.lat); ...
                  -(d0(dkdbt)+dPFDdbt)*DSIx(:,par.pindx.lat)];
        
        Sixx(:,kk) = mfactor(FD, tmp);
        kk = kk + 1;
    end
    % at aa
    if (par.opt_at == on & par.opt_aa == on)
        tmp = at*[d0(dkdat)*DSIx(:,par.pindx.aa); ...
                  -(d0(dkdat)+dPFDdat)*DSIx(:,par.pindx.aa)] + ...
              [-d0(Gx(:,par.pindx.lat))*dSi2Cdaa; ... %.*SIL; ...
               d0(Gx(:,par.pindx.lat))*dSi2Cdaa]; %.*SIL];
        
        Sixx(:,kk) = mfactor(FD, tmp);
        kk = kk + 1;
    end
    % at bb
    if (par.opt_at == on & par.opt_bb == on)
        tmp = at*[d0(dkdat)*DSIx(:,par.pindx.lbb); ...
                  -(d0(dkdat)+dPFDdat)*DSIx(:,par.pindx.lbb)] + ...
              [-d0(Gx(:,par.pindx.lat))*dSi2Cdbb; ... %.*SIL; ...
               d0(Gx(:,par.pindx.lat))*dSi2Cdbb]; %.*SIL];

        Sixx(:,kk) = mfactor(FD, tmp);
        kk = kk + 1;
    end
    % bt bt
    if (par.opt_bt == on)
        vout = buildPFD(par, parm);
        d2kdbtdbt = kappa_si./T.^2;
        d2PFDdbt2 = vout.d2PFDdbt2;
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
              [-d0(Gx(:,par.pindx.lbt))*dSi2Cdaa; ... %.*SIL; ...
               d0(Gx(:,par.pindx.lbt))*dSi2Cdaa]; %.*SIL];
        
        Sixx(:,kk) = mfactor(FD, tmp);
        kk = kk + 1;
    end
    % bt bb
    if (par.opt_bt == on & par.opt_bb == on)
        tmp = bt*[d0(dkdbt)*DSIx(:,par.pindx.lbb); ...
                  -(d0(dkdbt)+dPFDdbt)*DSIx(:,par.pindx.lbb)] + ...
              [-d0(Gx(:,par.pindx.lbt))*dSi2Cdbb; ... %.*SIL; ...
               d0(Gx(:,par.pindx.lbt))*dSi2Cdbb]; %.*SIL];

        Sixx(:,kk) = mfactor(FD, tmp);
        kk = kk + 1;
    end
    % aa aa
    if (par.opt_aa == on & par.opt_aa == on)
        tmp = [-Gxx(:,kk); ...
               Gxx(:,kk)];
        Sixx(:,kk) = mfactor(FD, tmp);
        kk = kk + 1;
    end
    % aa bb
    if (par.opt_aa == on & par.opt_bb == on)
        tmp = [-Gxx(:,kk); ...
               Gxx(:,kk)];
        Sixx(:,kk) = mfactor(FD, tmp);
        kk = kk + 1;
    end
    % bb bb
    if (par.opt_bb == on & par.opt_bb == on)
        tmp = [-Gxx(:,kk); ...
               Gxx(:,kk)];
        Sixx(:,kk) = mfactor(FD, tmp);
        kk = kk + 1;
    end
end
%
%% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function vout = buildPFD(par,parm);
M3d = parm.M3d;
grd = parm.grd;

[ny,nx,nz] = size(M3d);
M3D = zeros(ny,nx,nz+1);
M3D(:,:,1:end-1) = M3d;
% add the zw coordinate at the top of the extra layer
ZW3d = grd.ZW3d;
ZW3d = ZW3d(:,:,[1:end,end]);
ZW3d(:,:,end) = grd.ZW3d(:,:,end)+grd.dzt(end);
% areas of the top of the grid box
dAt = (grd.DXT3d.*grd.DYT3d).*M3d;
% volume of the grid boxes
dVt = (grd.DXT3d.*grd.DYT3d.*grd.DZT3d).*M3d;
%
n = nx*ny*(nz+1);
I0 = speye(n);
i0 = zeros(ny,nx,nz+1);
i0(:) = 1:n;
% periodic shifts OK because M3D has a layer of zeros on the bottom
iu = i0(:,:,[nz+1,1:nz]); %use a periodic upward shift
ib = i0(:,:,[2:nz+1,1]); % use a periodic downward shift
IU = I0(iu,:);
IB = I0(ib,:);
% keep only wet boxes
iocn = find(M3D(:));
I0 = I0(iocn,:); I0 = I0(:,iocn);
IU = IU(:,iocn); IU = IU(iocn,:);
IB = IB(:,iocn); IB = IB(iocn,:);
% (averages POP onto the top of the grid boxes)
AVG = d0((I0+IU)*M3D(iocn))\(I0+IU);
% (compute the divergence in the center of the grid boxes)
DIV = d0(dVt(iocn))\(I0-IB)*d0(dAt(iocn));
% (compute the flux at the top of the grid cells)
% mimics a Martin curve flux attenuation profile
%(see Kriest and Oschelies 2008 in Biogeosciences)
% ++++++++++++++++++++++++++++++++++++++++++++++++++
at  = parm.at;
bt  = parm.bt;
T  = parm.theta + 273.15;
kappa_si = at*exp(-bt./T); 
% ++++++++++++++++++++++++++++++++++++++++++++++++++

r = kappa_si;
b = parm.bsi;
a = r./b;
% particle sinking velocity at the top of the grid cells.
MSK = M3D.*M3D(:,:,[nz+1,1:nz]);
M = MSK.*ZW3d;

w = -a.*M(iocn);
%
dadb = -r./(b.^2);
dadr = 1./b;
drdat = exp(-bt./T);
drdbt = -(at*exp(-bt./T))./T;
dadat =  dadr.*drdat;
dadbt =  dadr.*drdbt;
%
dwdb = -dadb.*M(iocn);
dwdat = -dadat.*M(iocn);
dwdbt = -dadbt.*M(iocn);
%
d2adb2 = 2*r./(b^3);
d2adat2 = 0;
d2adbt2 = dadr.*(at*exp(-bt./T))./(T.^2);
d2adatdbt = -dadr.*exp(-bt./T)./T;
d2adatdb  = -(1./b.^2).*drdat;
d2adbtdb  = -(1./b.^2).*drdbt;

d2wdb2 = -d2adb2.*M(iocn);
d2wdat2 = 0;
d2wdbt2 = -d2adbt2.*M(iocn);
d2wdatdbt = -d2adatdbt.*M(iocn);
d2wdatdb  = -d2adatdb.*M(iocn);
d2wdbtdb  = -d2adbtdb.*M(iocn);
%FLUX = d0(w(iocn))*AVG;
FLUX = d0(w)*IU;
dFLUXdb = d0(dwdb)*IU;
dFLUXdat = d0(dwdat)*IU;
dFLUXdbt = d0(dwdbt)*IU;
%
d2FLUXdb2 = d0(d2wdb2)*IU;
d2FLUXdat2 = d0(d2wdat2)*IU;
d2FLUXdbt2 = d0(d2wdbt2)*IU;
d2FLUXdatdbt = d0(d2wdatdbt)*IU;
d2FLUXdatdb  = d0(d2wdatdb)*IU;
d2FLUXdbtdb  = d0(d2wdbtdb)*IU;
% particle flux divergence operator
vout.PFdiv = DIV*FLUX;
vout.dPFDdb = DIV*dFLUXdb;
vout.dPFDdat = DIV*dFLUXdat;
vout.dPFDdbt = DIV*dFLUXdbt;
%
vout.d2PFDdb2 = DIV*d2FLUXdb2;
vout.d2PFDdat2 = DIV*d2FLUXdat2;
vout.d2PFDdbt2 = DIV*d2FLUXdbt2;
vout.d2PFDdatdbt = DIV*d2FLUXdatdbt;
vout.d2PFDdatdb = DIV*d2FLUXdatdb;
vout.d2PFDdbtdb = DIV*d2FLUXdbtdb;
%% +++++++++++++++++++++++++++++++++++++++++++