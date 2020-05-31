function [f, fx] = neglogpost(x, parm, par)
on = true; off = false;
% reset parameters if they are too large/small;
[x] = reset_par(x, parm, par);
nx = length(x); % number of parameters

dVt  = parm.dVt  ;
M3d  = parm.M3d  ;
iwet = parm.iwet ;
nwet = parm.nwet ;

W   = d0(dVt(iwet)/sum(dVt(iwet)));
mu_dip = sum(W*parm.po4obs(iwet))/sum(diag(W));
var_dip = sum(W*(parm.po4obs(iwet)-mu_dip).^2)/sum(diag(W));
Wp = W/var_dip;
%
mu_dic = sum(W*parm.DICobs(iwet))/sum(diag(W));
var_dic = sum(W*(parm.DICobs(iwet)-mu_dic).^2)/sum(diag(W));
Wc = W/var_dic;
%
mu_o2 = sum(W*parm.o2obs)/sum(diag(W));
var_o2 = sum(W*(parm.o2obs-mu_o2).^2)/sum(diag(W));
Wo = W/var_o2;

%%%%%%%%%%%%%%%%%%   Slove P    %%%%%%%%%%%%%%%%%%%%%%%%
[P, Px, parm] = eqPcycle(par, parm, x);
DIP = M3d+nan;  DIP(iwet) = P(1+0*nwet:1*nwet) ;
POP = M3d+nan;  POP(iwet) = P(1+1*nwet:2*nwet) ;
DOP = M3d+nan;  DOP(iwet) = P(1+2*nwet:3*nwet) ;

parm.DIP = DIP(iwet); parm.Px = Px;
%%%%%%%%%%%%%%%%%%   End Solve P    %%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%%%%%%%%%%%   Slove C    %%%%%%%%%%%%%%%%%%%%%%%%
[C, Cx, parm] = eqCcycle(par, parm, x);
DIC = M3d+nan; DIC(iwet) = C(0*nwet+1:1*nwet) ;
POC = M3d+nan; POC(iwet) = C(1*nwet+1:2*nwet) ;
DOC = M3d+nan; DOC(iwet) = C(2*nwet+1:3*nwet) ;
CaC = M3d+nan; CaC(iwet) = C(3*nwet+1:4*nwet) ;
parm.DIC = DIC(iwet);      parm.DOC = DOC(iwet);
parm.DICx = Cx(1:nwet,:);  parm.DOCx = Cx(2*nwet+1:3*nwet,:);
%%%%%%%%%%%%%%%%%%   End Solve C    %%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%%%%%%%%%%%   Slove C    %%%%%%%%%%%%%%%%%%%%%%%%
[O, Ox] = eqOcycle(par, parm, x);
O2 = M3d+nan;  O2(iwet) = O;
%%%%%%%%%%%%%%%%%%   End Solve C    %%%%%%%%%%%%%%%%%%%%
%
DIC = DIC + parm.human_co2;
ep  = DIP(iwet) - parm.po4obs(iwet);
ec  = DIC(iwet) - parm.DICobs(iwet);
eo  = O - parm.o2obs;
%
f   = 0.5*(ep.'*Wp*ep + ec.'*Wc*ec + eo.'*Wo*eo);
fprintf('current objective function value is %3.3e \n',f);
%
if (nargout > 1)
    if (par.biogeochem.opt_sigma == on)
        px(:, par.pindx.lsigma) = Px(1:nwet,par.pindx.lsigma);
        cx(:, par.pindx.lsigma) = Cx(1:nwet,par.pindx.lsigma);
        ox(:, par.pindx.lsigma) = Ox(:, par.pindx.lsigma);
    end

    if (par.biogeochem.opt_slopep == on)
        px(:, par.pindx.slopep) = Px(1:nwet,par.pindx.slopep);
        cx(:, par.pindx.slopep) = Cx(1:nwet,par.pindx.slopep);
        ox(:, par.pindx.slopep) = Ox(:, par.pindx.slopep);
    end
    
    if (par.biogeochem.opt_interpp == on)
        px(:, par.pindx.linterpp) = Px(1:nwet,par.pindx.linterpp);
        cx(:, par.pindx.linterpp) = Cx(1:nwet,par.pindx.linterpp);
        ox(:, par.pindx.linterpp) = Ox(:, par.pindx.linterpp);
    end
    
    if (par.biogeochem.opt_alpha == on)
        px(:, par.pindx.lalpha) = Px(1:nwet,par.pindx.lalpha);
        cx(:, par.pindx.lalpha) = Cx(1:nwet,par.pindx.lalpha);
        ox(:, par.pindx.lalpha) = Ox(:, par.pindx.lalpha);
    end
    
    if (par.biogeochem.opt_beta == on)
        px(:, par.pindx.lbeta) = Px(1:nwet,par.pindx.lbeta);
        cx(:, par.pindx.lbeta) = Cx(1:nwet,par.pindx.lbeta);
        ox(:, par.pindx.lbeta) = Ox(:, par.pindx.lbeta);
    end
    
    if (par.biogeochem.opt_kappa_dp == on)
        px(:, par.pindx.lkappa_dp) = Px(1:nwet,par.pindx.lkappa_dp);
        cx(:, par.pindx.lkappa_dp) = Cx(1:nwet,par.pindx.lkappa_dp);
        ox(:, par.pindx.lkappa_dp) = Ox(:, par.pindx.lkappa_dp);
    end
    
    if (par.biogeochem.opt_slopec == on)
        px(:, par.pindx.slopec) = zeros(nwet,1);
        cx(:, par.pindx.slopec) = Cx(1:nwet,par.pindx.slopec);
        ox(:, par.pindx.slopec) = Ox(:, par.pindx.slopec);
    end
    
    if (par.biogeochem.opt_interpc == on)
        px(:, par.pindx.linterpc) = zeros(nwet,1);
        cx(:, par.pindx.linterpc) = Cx(1:nwet,par.pindx.linterpc);
        ox(:, par.pindx.linterpc) = Ox(:, par.pindx.linterpc);
    end
    
    if (par.biogeochem.opt_d == on)
        px(:, par.pindx.ld) = zeros(nwet,1);
        cx(:, par.pindx.ld) = Cx(1:nwet,par.pindx.ld);
        ox(:, par.pindx.ld) = Ox(:, par.pindx.ld);
    end
    
    if (par.biogeochem.opt_kappa_dc == on)
        px(:, par.pindx.lkappa_dc) = zeros(nwet,1);
        cx(:, par.pindx.lkappa_dc) = Cx(1:nwet,par.pindx.lkappa_dc);
        ox(:, par.pindx.lkappa_dc) = Ox(:, par.pindx.lkappa_dc);
    end
    
    if (par.biogeochem.opt_RR == on)
        px(:, par.pindx.lRR) = zeros(nwet,1);
        cx(:, par.pindx.lRR) = Cx(1:nwet,par.pindx.lRR);
        ox(:, par.pindx.lRR) = Ox(:, par.pindx.lRR);
    end

    if (par.biogeochem.opt_slopeo == on)
        px(:, par.pindx.slopeo) = zeros(nwet,1);
        cx(:, par.pindx.slopeo) = zeros(nwet,1);
        ox(:, par.pindx.slopeo) = Ox(:, par.pindx.slopeo);
    end

    if (par.biogeochem.opt_interpo == on)
        px(:, par.pindx.linterpo) = zeros(nwet,1);
        cx(:, par.pindx.linterpo) = zeros(nwet,1);
        ox(:, par.pindx.linterpo) = Ox(:,par.pindx.linterpo);
    end
    
    fx = ep.'*Wp*px + ec.'*Wc*cx + eo.'*Wo*ox;
end

fxx = zeros(4,4);
if (nargout>2)
    pxx = Pxx(1:nwet,:);
    % sigma sigma
    kk = 1;
    if (par.biogeochem.opt_sigma == on)
        fxx(par.pindx.lsigma, par.pindx.lsigma) = ...
            px(:,par.pindx.lsigma).'*W*px(:,par.pindx.lsigma) + ...
            ep.'*W*pxx(:,kk);
        kk = kk + 1;
    end

    % sigma kappa_dp
    if (par.biogeochem.opt_sigma == on & par.biogeochem.opt_kappa_dp == on)
        fxx(par.pindx.lsigma, par.pindx.lkappa_dp) = ...
            px(:,par.pindx.lsigma).'*W*px(:,par.pindx.lkappa_dp) + ...
            ep.'*W*pxx(:,kk);
        kk = kk + 1;
    end

    % sigma slopep
    if (par.biogeochem.opt_sigma == on & par.biogeochem.opt_slopep == on)
        fxx(par.pindx.lsigma, par.pindx.slopep) = ...
            px(:,par.pindx.lsigma).'*W*px(:,par.pindx.slopep) + ...
            ep.'*W*pxx(:,kk);
        kk = kk + 1;
    end

    % sigma interpp
    if (par.biogeochem.opt_sigma == on & par.biogeochem.opt_interpp == on)
        fxx(par.pindx.lsigma, par.pindx.linterpp) = ...
            px(:,par.pindx.lsigma).'*W*px(:,par.pindx.linterpp) + ...
            ep.'*W*pxx(:,kk);
        kk = kk + 1;
    end

    % sigma alpha
    if (par.biogeochem.opt_sigma == on & par.biogeochem.opt_alpha == on)
        fxx(par.pindx.lsigma, par.pindx.lalpha) = ...
            px(:,par.pindx.lsigma).'*W*px(:,par.pindx.lalpha) + ...
            ep.'*W*pxx(:,kk);
        kk = kk + 1;
    end

    % sigma beta
    if (par.biogeochem.opt_sigma == on & par.biogeochem.opt_beta == on)
        fxx(par.pindx.lsigma, par.pindx.lbeta) = ...
            px(:,par.pindx.lsigma).'*W*px(:,par.pindx.lbeta) + ...
            ep.'*W*pxx(:,kk);
        kk = kk + 1;
    end

    % kappa_dp kappa_dp
    if (par.biogeochem.opt_kappa_dp == on)
        fxx(par.pindx.lkappa_dp, par.pindx.lkappa_dp) = ...
            px(:,par.pindx.lkappa_dp).'*W*px(:,par.pindx.lkappa_dp) + ...
            ep.'*W*pxx(:,kk);
        kk = kk + 1;
    end

    % kappa_dp slopep
    if (par.biogeochem.opt_kappa_dp == on & par.biogeochem.opt_slopep == on)
        fxx(par.pindx.lkappa_dp, par.pindx.slopep) = ...
            px(:,par.pindx.lkappa_dp).'*W*px(:,par.pindx.slopep) + ...
            ep.'*W*pxx(:,kk);
        kk = kk + 1;
    end

    % kappa_dp interpp
    if (par.biogeochem.opt_kappa_dp == on & par.biogeochem.opt_interpp == on)
        fxx(par.pindx.lkappa_dp, par.pindx.linterpp) = ...
            px(:,par.pindx.lkappa_dp).'*W*px(:,par.pindx.linterpp) + ...
            ep.'*W*pxx(:,kk);
        kk = kk + 1;
    end

    % kappa_dp alpha
    if (par.biogeochem.opt_kappa_dp == on & par.biogeochem.opt_alpha == on)
        fxx(par.pindx.lkappa_dp, par.pindx.lalpha) = ...
            px(:,par.pindx.lkappa_dp).'*W*px(:,par.pindx.lalpha) + ...
            ep.'*W*pxx(:,kk);
        kk = kk + 1;
    end

    % kappa_dp beta
    if (par.biogeochem.opt_kappa_dp == on & par.biogeochem.opt_beta == on)
        fxx(par.pindx.lkappa_dp, par.pindx.lbeta) = ...
            px(:,par.pindx.lkappa_dp).'*W*px(:,par.pindx.lbeta) + ...
            ep.'*W*pxx(:,kk);
        kk = kk + 1;
    end

    % slopep slopep
    if (par.biogeochem.opt_slopep == on)
        fxx(par.pindx.slopep, par.pindx.slopep) = ...
            px(:,par.pindx.slopep).'*W*px(:,par.pindx.slopep) + ...
            ep.'*W*pxx(:,kk);
        kk = kk + 1;
    end
    
    % slopep interpp
    if (par.biogeochem.opt_slopep == on & par.biogeochem.opt_interpp ...
        == on)
        fxx(par.pindx.slopep, par.pindx.linterpp) = ...
            px(:,par.pindx.slopep).'*W*px(:,par.pindx.linterpp) + ...
            ep.'*W*pxx(:,kk);
        kk = kk + 1;
    end

    % slopep alpha
    if (par.biogeochem.opt_slopep == on & par.biogeochem.opt_alpha == on)
        fxx(par.pindx.slopep, par.pindx.lalpha) = ...
            px(:,par.pindx.slopep).'*W*px(:,par.pindx.lalpha) + ...
            ep.'*W*pxx(:,kk);
        kk = kk + 1;
    end

    % slopep beta
    if (par.biogeochem.opt_slopep == on & par.biogeochem.opt_beta == on)
        fxx(par.pindx.slopep, par.pindx.lbeta) = ...
            px(:,par.pindx.slopep).'*W*px(:,par.pindx.lbeta) + ...
            ep.'*W*pxx(:,kk);
        kk = kk + 1;
    end
    
    % interpp interpp
    if (par.biogeochem.opt_interpp == on)
        fxx(par.pindx.linterpp, par.pindx.linterpp) = ...
            px(:,par.pindx.linterpp).'*W*px(:,par.pindx.linterpp) + ...
            ep.'*W*pxx(:,kk);
        kk = kk + 1;
    end    

    % interpp alpha
    if (par.biogeochem.opt_interpp == on & par.biogeochem.opt_alpha == on)
        fxx(par.pindx.linterpp, par.pindx.lalpha) = ...
            px(:,par.pindx.linterpp).'*W*px(:,par.pindx.lalpha) + ...
            ep.'*W*pxx(:,kk);
        kk = kk + 1;
    end

    % interpp beta
    if (par.biogeochem.opt_interpp == on & par.biogeochem.opt_beta == on)
        fxx(par.pindx.linterpp, par.pindx.lbeta) = ...
            px(:,par.pindx.linterpp).'*W*px(:,par.pindx.lbeta) + ...
            ep.'*W*pxx(:,kk);
        kk = kk + 1;
    end
    
    % alpha alpha
    if (par.biogeochem.opt_alpha == on)
        fxx(par.pindx.lalpha, par.pindx.lalpha) = ...
            px(:,par.pindx.lalpha).'*W*px(:,par.pindx.lalpha) + ...
            ep.'*W*pxx(:,kk);
        kk = kk + 1;
    end

    % alpha beta
    if (par.biogeochem.opt_alpha == on & par.biogeochem.opt_beta == on)
        fxx(par.pindx.lalpha, par.pindx.lbeta) = ...
            px(:,par.pindx.lalpha).'*W*px(:,par.pindx.lbeta) + ...
            ep.'*W*pxx(:,kk);
        kk = kk + 1;
    end
    
    % beta beta
    if (par.biogeochem.opt_beta == on)
        fxx(par.pindx.lbeta, par.pindx.lbeta) = ...
            px(:,par.pindx.lbeta).'*W*px(:,par.pindx.lbeta) + ...
            ep.'*W*pxx(:,kk);
        kk = kk + 1;
    end
end

fxx = 0.5*(fxx + fxx.'); % symetric