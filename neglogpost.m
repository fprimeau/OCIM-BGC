function [f, fx, fxx] = neglogpost(x, parm, par)
on = true; off = false;
% reset parameters if they are too large/small;
[x] = reset_par(x, parm, par);
nx = length(x); % number of parameters

dVt  = parm.dVt  ;
M3d  = parm.M3d  ;
iwet = parm.iwet ;
nwet = parm.nwet ;
%
f = 0;
%%%%%%%%%%%%%%%%%%   Slove P    %%%%%%%%%%%%%%%%%%%%%%%%
%
W   = d0(dVt(iwet)/sum(dVt(iwet)));
mu_dip  = sum(W*parm.po4obs(iwet))/sum(diag(W));
var_dip = sum(W*(parm.po4obs(iwet)-mu_dip).^2)/sum(diag(W));
Wp = W/var_dip;
%
[parm, P, Px, Pxx] = eqPcycle(par, parm, x);
DIP = M3d+nan;  DIP(iwet) = P(1+0*nwet:1*nwet) ;
POP = M3d+nan;  POP(iwet) = P(1+1*nwet:2*nwet) ;
DOP = M3d+nan;  DOP(iwet) = P(1+2*nwet:3*nwet) ;
parm.DIP = DIP(iwet);  parm.Px = Px;
% parm.Pxx = Pxx;
% DIP error
ep = DIP(iwet) - parm.po4obs(iwet);
f = f + 0.5*(ep.'*Wp*ep); 
%%%%%%%%%%%%%%%%%%   End Solve P    %%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%%%%%%%%%%%   End Solve Si    %%%%%%%%%%%%%%%%%%%%
if (par.Simodel == on)
    %
    mu_sil = sum(W*parm.SIL(iwet))/sum(diag(W));
    var_sil = sum(W*(parm.SIL(iwet)-mu_sil).^2)/sum(diag(W));
    Ws = W/var_sil;
    %
    [Si,Six] = eqSicycle(par, parm, x);
    SIO = M3d+nan;  SIO(iwet) = Si(1:nwet);
    DSI = M3d+nan;  DSI(iwet) = Si(nwet+1:end);
    % SiO error
    es = SIO(iwet) - parm.SIL(iwet);
    f = f + 0.5*(es.'*Ws*es); 
end
%%%%%%%%%%%%%%%%%%   End Solve Si    %%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%%%%%%%%%%%   Slove C    %%%%%%%%%%%%%%%%%%%%%%%%
if (par.Cmodel == on)
    mu_dic = sum(W*parm.DICobs(iwet))/sum(diag(W));
    var_dic = sum(W*(parm.DICobs(iwet)-mu_dic).^2)/sum(diag(W));
    Wc = W/var_dic;
    %
    [C, Cx, parm] = eqCcycle(par, parm, x);
    DIC = M3d+nan; DIC(iwet) = C(0*nwet+1:1*nwet) ;
    POC = M3d+nan; POC(iwet) = C(1*nwet+1:2*nwet) ;
    DOC = M3d+nan; DOC(iwet) = C(2*nwet+1:3*nwet) ;
    CaC = M3d+nan; CaC(iwet) = C(3*nwet+1:4*nwet) ;
    parm.DIC = DIC(iwet);      parm.DOC = DOC(iwet);
    parm.DICx = Cx(1:nwet,:);  parm.DOCx = Cx(2*nwet+1:3*nwet,:);
    % DIC error
    DIC = DIC + parm.human_co2;
    ec  = DIC(iwet) - parm.DICobs(iwet);
    f   = f + 0.5*(ec.'*Wc*ec);
end
%%%%%%%%%%%%%%%%%%   End Solve C    %%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%%%%%%%%%%%   Slove O    %%%%%%%%%%%%%%%%%%%%%%%%
if (par.Omodel == on)
    %
    mu_o2 = sum(W*parm.o2obs)/sum(diag(W));
    var_o2 = sum(W*(parm.o2obs-mu_o2).^2)/sum(diag(W));
    Wo = W/var_o2;
    %
    [O, Ox] = eqOcycle(par, parm, x);
    O2 = M3d+nan;  O2(iwet) = O;
    eo  = O - parm.o2obs;
    f   = f + 0.5*(eo.'*Wo*eo);
end
%%%%%%%%%%%%%%%%%%   End Solve O    %%%%%%%%%%%%%%%%%%%%
fprintf('current objective function value is %3.3e \n',f);
%
%% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
if (nargout > 1)
    fx = zeros(length(x), 1);
    if (par.opt_sigma == on)
        fx(par.pindx.lsigma) = ep.'*Wp*Px(1:nwet,par.pindx.lsigma);
        if (par.Cmodel == on)
            fx(par.pindx.lsigma) = fx(par.pindx.lsigma) + ...
                ec.'*Wc*Cx(1:nwet,par.pindx.lsigma);
        end 
        if (par.Simodel) == on
            fx(par.pindx.lsigma) = fx(par.pindx.lsigma) + ...
                es.'*Ws*Six(1:nwet,par.pindx.lsigma);
        end 
        if (par.Omodel) == on
            ox(:, par.pindx.lsigma) = fx(par.pindx.lsigma) + ...
                eo.'*Wo*Ox(:, par.pindx.lsigma);
        end
    end

    if (par.opt_slopep == on)
        fx(par.pindx.slopep) = ep.'*Wp*Px(1:nwet,par.pindx.slopep);
        if (par.Cmodel == on)
            fx(par.pindx.slopep) = fx(par.pindx.slopep) + ...
                ec.'*Wc*Cx(1:nwet,par.pindx.slopep);
        end 
        if (par.Simodel == on) 
            fx(par.pindx.slopep) = fx(par.pindx.slopep) + ...
                es.'*Ws*Six(1:nwet,par.pindx.slopep);
        end 
        if (par.Omodel == on) 
            fx(par.pindx.slopep) = fx(par.pindx.slopep) + ...
                eo.'*Wo*Ox(:, par.pindx.slopep);
        end 
    end
    
    if (par.opt_interpp == on)
        fx(par.pindx.linterpp) = ep.'*Wp*Px(1:nwet,par.pindx.linterpp);
        if (par.Cmodel == on) 
            cx(:, par.pindx.linterpp) = fx(par.pindx.linterpp) + ...
                ec.'*Wc*Cx(1:nwet,par.pindx.linterpp);
        end 
        if (par.Simodel == on)
            fx(par.pindx.linterpp) = fx(par.pindx.linterpp) + ...
                es.'*Ws*Six(1:nwet, par.pindx.linterpp);
        end 
        if (par.Omodel == on) 
            fx(par.pindx.linterpp) = fx(par.pindx.linterpp) + ...
                eo.'*Wo*Ox(:, par.pindx.linterpp);
        end 
    end
    
    if (par.opt_alpha == on)
        fx(par.pindx.lalpha) = ep.'*Wp*Px(1:nwet,par.pindx.lalpha);
        if (par.Cmodel == on) 
            fx(par.pindx.lalpha) = fx(par.pindx.lalpha) + ...
                ec.'*Wc*Cx(1:nwet,par.pindx.lalpha);
        end 
        if (par.Simodel == on) 
            fx(par.pindx.lalpha) = fx(par.pindx.lalpha) + ...
                es.'*Ws*Six(1:nwet,par.pindx.lalpha);
        end 
        if (par.Omodel == on)
            fx(par.pindx.lalpha) = fx(par.pindx.lalpha) + ...
                eo.'*Wo*Ox(:, par.pindx.lalpha);
        end 
    end
    
    if (par.opt_beta == on)
        fx(par.pindx.lbeta) = ep.'*Wp*Px(1:nwet,par.pindx.lbeta);
        if (par.Cmodel == on) 
            fx(par.pindx.lbeta) = fx(par.pindx.lbeta) + ...
                ec.'*Wc*Cx(1:nwet,par.pindx.lbeta);
        end 
        if (par.Simodel == on)
            fx(par.pindx.lbeta) = fx(par.pindx.lbeta) + ...
                es.'*Ws*Six(1:nwet,par.pindx.lbeta);
        end 
        if (par.Omodel == on) 
            fx(par.pindx.lbeta) = fx(par.pindx.lbeta) + ...
                eo.'*Wo*Ox(:, par.pindx.lbeta);
        end 
    end
    
    if (par.opt_kappa_dp == on)
        fx(par.pindx.lkappa_dp) = ep.'*Wp*Px(1:nwet,par.pindx.lkappa_dp);
        if (par.Cmodel == on) 
            fx(par.pindx.lkappa_dp) = fx(par.pindx.lkappa_dp) + ...
                ec.'*Wc*Cx(1:nwet,par.pindx.lkappa_dp);
        end
        if (par.Simodel == on) 
            fx(par.pindx.lkappa_dp) = fx(par.pindx.lkappa_dp) + ...
                es.'*Ws*Six(1:nwet,par.pindx.lkappa_dp);
        end
        if (par.Omodel == on) 
            fx(par.pindx.lkappa_dp) = fx(par.pindx.lkappa_dp) + ...
                eo.'*Wo*Ox(:, par.pindx.lkappa_dp);
        end 
    end
    
    if (par.opt_slopec == on)
        if (par.Cmodel == on) 
            fx(par.pindx.slopec) = fx(par.pindx.slopec) + ...
                ec.'*Wc*Cx(1:nwet,par.pindx.slopec);
        end 
        if (par.Omodel == on) 
            fx(par.pindx.slopec) = fx(par.pindx.slopec) + ...
                eo.'*Wo*Ox(:, par.pindx.slopec);
        end 
    end
    
    if (par.opt_interpc == on)
        if (par.Cmodel == on) 
            fx(par.pindx.linterpc) = ec.'*Wc*Cx(1:nwet, ...
                                                par.pindx.linterpc);
        end 
        if (par.Omodel == on)
            fx(par.pindx.linterpc) = fx(par.pindx.linterpc) + ...
                eo.'*Wo*Ox(:, par.pindx.linterpc);
        end 
    end
    
    if (par.opt_d == on)
        if (par.Cmodel == on) 
            fx(par.pindx.ld) = ec.'*Wc*Cx(1:nwet,par.pindx.ld);
        end 
        if (par.Omodel == on) 
            fx(par.pindx.ld) = fx(par.pindx.ld) + ...
                eo.'*Wo*Ox(:, par.pindx.ld);
        end 
    end
    
    if (par.opt_kappa_dc == on)
        if (par.Cmodel == on) 
            fx(par.pindx.lkappa_dc) = ec.'*Wc*Cx(1:nwet, ...
                                                 par.pindx.lkappa_dc);
        end 
        if (par.Omodel == on) 
            fx(par.pindx.lkappa_dc) = fx(par.pindx.lkappa_dc) + ...
                eo.'*Wo*Ox(:, par.pindx.lkappa_dc);
        end 
    end
    
    if (par.opt_RR == on)
        if (par.Cmodel == on)
            fx(par.pindx.lRR) = ec.'*Wc*Cx(1:nwet, par.pindx.lRR);
        end 
        if (par.Omodel == on) 
            fx(par.pindx.lRR) = fx(par.pindx.lRR) + ...
                eo.'*Wo*Ox(:, par.pindx.lRR);
        end 
    end

    if (par.Omodel == on)
        if (par.opt_slopeo == on)
            fx(par.pindx.slopeo) = eo.'*Wo*Ox(:, par.pindx.slopeo);
        end 
    end

    if (par.Omodel == on)
        if (par.opt_interpo == on)
            fx(par.pindx.linterpo) = eo.'*Wo*Ox(:, par.pindx.linterpo);
        end 
    end

    if (par.Simodel == on)
        if (par.opt_at == on)
            fx(par.pindx.lat) = es.'*Ws*Six(1:nwet, par.pindx.lat);
        end 
    end

    if (par.Simodel == on) 
        if (par.opt_bt == on)
            fx(par.pindx.lbt) = es.'*Ws*Six(1:nwet, par.pindx.lbt);
        end 
    end

    if (par.Simodel == on) 
            if (par.opt_aa == on)
            fx(par.pindx.aa) = es.'*Ws*Six(1:nwet, par.pindx.aa);
        end  
    end

    if (par.Simodel == on)
        if (par.opt_bb == on)
        fx(par.pindx.lbb) = es.'*Ws*Six(1:nwet, par.pindx.lbb);
        end 
    end

    if (par.Simodel == on)
        if (par.opt_kappa_gs == on)
            fx(par.pindx.lkappa_gs) = es.'*Ws*Six(1:nwet, ...
                                                  par.pindx.lkappa_gs);
        end  
    end
end

%% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
if (nargout>2)
    fxx = zeros(4,4);
    px  = Px(1:nwet,:);
    pxx = Pxx(1:nwet,:);
    % sigma sigma
    kk = 1;
    if (par.opt_sigma == on)
        fxx(par.pindx.lsigma, par.pindx.lsigma) = ...
            px(:,par.pindx.lsigma).'*Wp*px(:,par.pindx.lsigma) + ...
            ep.'*Wp*pxx(:,kk);
        fxx(par.pindx.lsigma, par.pindx.lsigma) = fxx(par.pindx.lsigma, ...
                                                      par.pindx.lsigma);
        kk = kk + 1;
    end

    % sigma kappa_dp
    if (par.opt_sigma == on & par.opt_kappa_dp == on)
        fxx(par.pindx.lsigma, par.pindx.lkappa_dp) = ...
            px(:,par.pindx.lsigma).'*Wp*px(:,par.pindx.lkappa_dp) + ...
            ep.'*Wp*pxx(:,kk);
        fxx(par.pindx.lkappa_dp, par.pindx.lsigma) = fxx(par.pindx.lsigma, ...
                                                         par.pindx.lkappa_dp);
        kk = kk + 1;
    end

    % sigma slopep
    if (par.opt_sigma == on & par.opt_slopep == on)
        fxx(par.pindx.lsigma, par.pindx.slopep) = ...
            px(:,par.pindx.lsigma).'*Wp*px(:,par.pindx.slopep) + ...
            ep.'*Wp*pxx(:,kk);
        fxx(par.pindx.slopep, par.pindx.lsigma) = fxx(par.pindx.lsigma, ...
                                                      par.pindx.slopep);
        
        kk = kk + 1;
    end

    % sigma interpp
    if (par.opt_sigma == on & par.opt_interpp == on)
        fxx(par.pindx.lsigma, par.pindx.linterpp) = ...
            px(:,par.pindx.lsigma).'*Wp*px(:,par.pindx.linterpp) + ...
            ep.'*Wp*pxx(:,kk);
        fxx(par.pindx.linterpp,par.pindx.lsigma) = fxx(par.pindx.lsigma, ...
                                                       par.pindx.linterpp);
        kk = kk + 1;
    end

    % sigma alpha
    if (par.opt_sigma == on & par.opt_alpha == on)
        fxx(par.pindx.lsigma, par.pindx.lalpha) = ...
            px(:,par.pindx.lsigma).'*Wp*px(:,par.pindx.lalpha) + ...
            ep.'*Wp*pxx(:,kk);
        fxx(par.pindx.lalpha, par.pindx.lsigma) = fxx(par.pindx.lsigma, ...
                                                      par.pindx.lalpha);
        kk = kk + 1;
    end

    % sigma beta
    if (par.opt_sigma == on & par.opt_beta == on)
        fxx(par.pindx.lsigma, par.pindx.lbeta) = ...
            px(:,par.pindx.lsigma).'*Wp*px(:,par.pindx.lbeta) + ...
            ep.'*Wp*pxx(:,kk);
        fxx(par.pindx.lbeta, ...
            par.pindx.lsigma)= fxx(par.pindx.lsigma, par.pindx.lbeta);
        kk = kk + 1;
    end

    % kappa_dp kappa_dp
    if (par.opt_kappa_dp == on)
        fxx(par.pindx.lkappa_dp, par.pindx.lkappa_dp) = ...
            px(:,par.pindx.lkappa_dp).'*Wp*px(:,par.pindx.lkappa_dp) + ...
            ep.'*Wp*pxx(:,kk);
        fxx(par.pindx.lkappa_dp, par.pindx.lkappa_dp) = ...
            fxx(par.pindx.lkappa_dp, par.pindx.lkappa_dp);
        kk = kk + 1;
    end

    % kappa_dp slopep
    if (par.opt_kappa_dp == on & par.opt_slopep == on)
        fxx(par.pindx.lkappa_dp, par.pindx.slopep) = ...
            px(:,par.pindx.lkappa_dp).'*Wp*px(:,par.pindx.slopep) + ...
            ep.'*Wp*pxx(:,kk);
        fxx(par.pindx.slopep, par.pindx.lkappa_dp) = ...
            fxx(par.pindx.lkappa_dp, par.pindx.slopep);
        kk = kk + 1;
    end

    % kappa_dp interpp
    if (par.opt_kappa_dp == on & par.opt_interpp == on)
        fxx(par.pindx.lkappa_dp, par.pindx.linterpp) = ...
            px(:,par.pindx.lkappa_dp).'*Wp*px(:,par.pindx.linterpp) + ...
            ep.'*Wp*pxx(:,kk);
        fxx(par.pindx.linterpp, par.pindx.lkappa_dp) = ...
            fxx(par.pindx.lkappa_dp, par.pindx.linterpp);
        kk = kk + 1;
    end

    % kappa_dp alpha
    if (par.opt_kappa_dp == on & par.opt_alpha == on)
        fxx(par.pindx.lkappa_dp, par.pindx.lalpha) = ...
            px(:,par.pindx.lkappa_dp).'*Wp*px(:,par.pindx.lalpha) + ...
            ep.'*Wp*pxx(:,kk);
        fxx(par.pindx.lalpha, par.pindx.lkappa_dp) = ...
            fxx(par.pindx.lkappa_dp, par.pindx.lalpha);
        kk = kk + 1;
    end

    % kappa_dp beta
    if (par.opt_kappa_dp == on & par.opt_beta == on)
        fxx(par.pindx.lkappa_dp, par.pindx.lbeta) = ...
            px(:,par.pindx.lkappa_dp).'*Wp*px(:,par.pindx.lbeta) + ...
            ep.'*Wp*pxx(:,kk);
        fxx(par.pindx.lbeta, par.pindx.lkappa_dp) = ...
            fxx(par.pindx.lkappa_dp, par.pindx.lbeta);
        kk = kk + 1;
    end

    % slopep slopep
    if (par.opt_slopep == on)
        fxx(par.pindx.slopep, par.pindx.slopep) = ...
            px(:,par.pindx.slopep).'*Wp*px(:,par.pindx.slopep) + ...
            ep.'*Wp*pxx(:,kk);
        kk = kk + 1;
    end
    
    % slopep interpp
    if (par.opt_slopep == on & par.opt_interpp ...
        == on)
        fxx(par.pindx.slopep, par.pindx.linterpp) = ...
            px(:,par.pindx.slopep).'*Wp*px(:,par.pindx.linterpp) + ...
            ep.'*Wp*pxx(:,kk);
        fxx(par.pindx.linterpp, ...
            par.pindx.slopep)= fxx(par.pindx.slopep, par.pindx.linterpp);
        kk = kk + 1;
    end

    % slopep alpha
    if (par.opt_slopep == on & par.opt_alpha == on)
        fxx(par.pindx.slopep, par.pindx.lalpha) = ...
            px(:,par.pindx.slopep).'*Wp*px(:,par.pindx.lalpha) + ...
            ep.'*Wp*pxx(:,kk);
        fxx(par.pindx.lalpha, ...
            par.pindx.slopep)= fxx(par.pindx.slopep, par.pindx.lalpha);
        kk = kk + 1;
    end

    % slopep beta
    if (par.opt_slopep == on & par.opt_beta == on)
        fxx(par.pindx.slopep, par.pindx.lbeta) = ...
            px(:,par.pindx.slopep).'*Wp*px(:,par.pindx.lbeta) + ...
            ep.'*Wp*pxx(:,kk);
        fxx(par.pindx.lbeta, par.pindx.slopep) = fxx(par.pindx.slopep, ...
                                                     par.pindx.lbeta);
        kk = kk + 1;
    end
    
    % interpp interpp
    if (par.opt_interpp == on)
        fxx(par.pindx.linterpp, par.pindx.linterpp) = ...
            px(:,par.pindx.linterpp).'*Wp*px(:,par.pindx.linterpp) + ...
            ep.'*Wp*pxx(:,kk);
        fxx(par.pindx.linterpp, par.pindx.linterpp) = ...
            fxx(par.pindx.linterpp, par.pindx.linterpp);
        kk = kk + 1;
    end    

    % interpp alpha
    if (par.opt_interpp == on & par.opt_alpha == on)
        fxx(par.pindx.linterpp, par.pindx.lalpha) = ...
            px(:,par.pindx.linterpp).'*Wp*px(:,par.pindx.lalpha) + ...
            ep.'*Wp*pxx(:,kk);
        fxx(par.pindx.lalpha, ...
            par.pindx.linterpp)= fxx(par.pindx.linterpp, par.pindx.lalpha);
        kk = kk + 1;
    end

    % interpp beta
    if (par.opt_interpp == on & par.opt_beta == on)
        fxx(par.pindx.linterpp, par.pindx.lbeta) = ...
            px(:,par.pindx.linterpp).'*Wp*px(:,par.pindx.lbeta) + ...
            ep.'*Wp*pxx(:,kk);
        fxx(par.pindx.lbeta, ...
            par.pindx.linterpp)= fxx(par.pindx.linterpp, par.pindx.lbeta);
        kk = kk + 1;
    end
    
    % alpha alpha
    if (par.opt_alpha == on)
        fxx(par.pindx.lalpha, par.pindx.lalpha) = ...
            px(:,par.pindx.lalpha).'*Wp*px(:,par.pindx.lalpha) + ...
            ep.'*Wp*pxx(:,kk);
        kk = kk + 1;
    end

    % alpha beta
    if (par.opt_alpha == on & par.opt_beta == on)
        fxx(par.pindx.lalpha, par.pindx.lbeta) = ...
            px(:,par.pindx.lalpha).'*Wp*px(:,par.pindx.lbeta) + ...
            ep.'*Wp*pxx(:,kk);
        fxx(par.pindx.lbeta, ...
            par.pindx.lalpha)= fxx(par.pindx.lalpha, par.pindx.lbeta);
        kk = kk + 1;
    end
    
    % beta beta
    if (par.opt_beta == on)
        fxx(par.pindx.lbeta, par.pindx.lbeta) = ...
            px(:,par.pindx.lbeta).'*Wp*px(:,par.pindx.lbeta) + ...
            ep.'*Wp*pxx(:,kk);
        kk = kk + 1;
    end
    %% ---------------------------------------------------------------
    if par.Simodel == on
        if (par.opt_sigma == on)
            Sxx(:,kk) = sparse(2*nwet,1);
            kk = kk + 1;
        end

        % sigma kappa_dp
        if (par.opt_sigma == on & par.opt_kappa_dp == on)
            Sxx(:,kk) = sparse(2*nwet,1);
            kk = kk + 1;
        end

        % sigma slopep
        if (par.opt_sigma == on & par.opt_slopep == on)
            Sxx(:,kk) = sparse(2*nwet,1);
            kk = kk + 1;
        end

        % sigma interpp
        if (par.opt_sigma == on & par.opt_interpp == on)
            Sxx(:,kk) = sparse(2*nwet,1);
            kk = kk + 1;
        end

        % sigma alpha
        if (par.opt_sigma == on & par.opt_alpha == on)
            Sxx(:,kk) = sparse(2*nwet,1);
            kk = kk + 1;
        end

        % sigma beta
        if (par.opt_sigma == on & par.opt_beta == on)
            Sxx(:,kk) = sparse(2*nwet,1);
            kk = kk + 1;
        end

        % kappa_dp kappa_dp
        if (par.opt_kappa_dp == on)
            Sxx(:,kk) = sparse(2*nwet,1);
            kk = kk + 1;
        end

        % kappa_dp slopep
        if (par.opt_kappa_dp == on & par.opt_slopep == on)
            Sxx(:,kk) = sparse(2*nwet,1);
            kk = kk + 1;
        end

        % kappa_dp interpp
        if (par.opt_kappa_dp == on & par.opt_interpp == on)
            Sxx(:,kk) = sparse(2*nwet,1);
            kk = kk + 1;
        end

        % kappa_dp alpha
        if (par.opt_kappa_dp == on & par.opt_alpha == on)
            Sxx(:,kk) = sparse(2*nwet,1);
            kk = kk + 1;
        end

        % kappa_dp beta
        if (par.opt_kappa_dp == on & par.opt_beta == on)
            Sxx(:,kk) = sparse(2*nwet,1);
            kk = kk + 1;
        end

        % slopep slopep
        if (par.opt_slopep == on)
            Sxx(:,kk) = sparse(2*nwet,1);
            kk = kk + 1;
        end
        
        % slopep interpp
        if (par.opt_slopep == on & par.opt_interpp == on)
            Sxx(:,kk) = sparse(2*nwet,1);
            kk = kk + 1;
        end

        % slopep alpha
        if (par.opt_slopep == on & par.opt_alpha == on)
            Sxx(:,kk) = sparse(2*nwet,1);
            kk = kk + 1;
        end

        % slopep beta
        if (par.opt_slopep == on & par.opt_beta == on)
            Sxx(:,kk) = sparse(2*nwet,1);
            kk = kk + 1;
        end
        
        % interpp interpp
        if (par.opt_interpp == on)
            Sxx(:,kk) = sparse(2*nwet,1);
            kk = kk + 1;
        end    

        % interpp alpha
        if (par.opt_interpp == on & par.opt_alpha == on)
            Sxx(:,kk) = sparse(2*nwet,1);
            kk = kk + 1;
        end

        % interpp beta
        if (par.opt_interpp == on & par.opt_beta == on)
            Sxx(:,kk) = sparse(2*nwet,1);
            kk = kk + 1;
        end
        
        % alpha alpha
        if (par.opt_alpha == on)
            tmp = [-Si2C.*Gxx(:,kk); ...
                   Si2C.*Gxx(:,kk)]; %drhsdalphadalpha
            
            Sxx(:,kk) = mfactor(FFp, tmp);
            kk = kk + 1;
        end

        % alpha beta
        if (par.opt_alpha == on & par.opt_beta == on)
            tmp = [-Si2C.*Gxx(:,kk); ...
                   Si2C.*Gxx(:,kk)]; %drhsdalphadbeta
            
            Sxx(:,kk) = mfactor(FFp,-tmp);
            kk = kk + 1;
        end
        
        % beta beta
        if (par.opt_beta == on)
            tmp = [-Si2C.*Gxx(:,kk); ...
                   Si2C.*Gxx(:,kk)]; %drhsdalphadbeta
            
            Sxx(:,kk) = mfactor(FFp,-tmp);
            kk = kk + 1;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Si parameters
        % sigma at
        if (par.opt_sigma == on & par.opt_at == on)
            tmp = at*[-dkdat.*DSI(:,par.pindx.lsigma); ...
                      (d0(dkdat) + dPFDdat)*DSI(:,par.pindx.lsigma)];
            Sxx(:,kk) = mfactor(FFp,-tmp);
            kk = kk + 1;
        end
        % sigma bt
        if (par.opt_sigma == on & par.opt_bt == on)
            tmp = -bt*[-dkdbt.*DSI(:,par.pindx.lsigma); ...
                       (d0(dkdbt) + dPFDdbt)*DSI(:,par.pindx.lsigma)];
            Sxx(:,kk) = mfactor(FFp,-tmp);
            kk = kk + 1;
        end
        % sigma aa
        if (par.opt_sigma == on & par.opt_aa == on)
            tmp = -[d0(Gx(:,par.pindx.lsigma))*dSi2Cdaa; ... %.*SIL; ...
                    -d0(G(:,par.pindx.lsigma))*dSi2Cdaa]; %.*SIL];
            
            Sxx(:,kk) = mfactor(FD, tmp);
        end
        % sigma bb
        if (par.opt_sigma == on & par.opt_bb == on)
            tmp = -bb*[d0(Gx(par.pindx.lbb))*dSi2Cdbb; ... %.*SIL; ...
                       -d0(Gx(par.pindx.lbb))*dSi2Cdbb]; %.*SIL];

            Sxx(:,kk) = mfactor(FD, tmp);
        end
        % kappa_dp at
        if (par.opt_kappa_dp == on & par.opt_at == on)
            tmp = at*[-dkdat.*DSI(:,par.pindx.lkappa_dp); ...
                      (d0(dkdat) + dPFDdat)*DSI(:,par.pindx.lkappa_dp)];
            Sxx(:,kk) = mfactor(FFp,-tmp);
            kk = kk + 1;
        end
        % kappa_dp bt
        if (par.opt_kappa_dp == on & par.opt_bt == on)
            tmp = -bt*[-dkdbt.*DSI(:,par.pindx.lkappa_dp); ...
                       (d0(dkdbt) + dPFDdbt)*DSI(:,par.pindx.lkappa_dp)];
            Sxx(:,kk) = mfactor(FFp,-tmp);
            kk = kk + 1;
        end
        % kappa_dp aa
        if (par.opt_kappa_dp == on & par.opt_aa == on)
            tmp = -[d0(Gx(:,par.pindx.lkappa_dp))*dSi2Cdaa; ... %.*SIL; ...
                    -d0(G(:,par.pindx.lkappa_dp))*dSi2Cdaa]; %.*SIL];
            
            Sxx(:,kk) = mfactor(FD, tmp);
        end
        % kappa_dp bb
        if (par.opt_kappa_dp == on & par.opt_bb == on)
            tmp = -bb*[d0(Gx(par.pindx.lbb))*dSi2Cdbb; ... %.*SIL; ...
                       -d0(Gx(par.pindx.lbb))*dSi2Cdbb]; %.*SIL];

            Sxx(:,kk) = mfactor(FD, tmp);
        end
        % slope at
        if (par.opt_slope == on & par.opt_at == on)
            tmp = at*[-dkdat.*DSI(:,par.pindx.slope); ...
                      (d0(dkdat) + dPFDdat)*DSI(:,par.pindx.slope)];
            Sxx(:,kk) = mfactor(FFp,-tmp);
            kk = kk + 1;
        end
        % slope bt
        if (par.opt_slope == on & par.opt_bt == on)
            tmp = -bt*[-dkdbt.*DSI(:,par.pindx.slope); ...
                       (d0(dkdbt) + dPFDdbt)*DSI(:,par.pindx.slope)];
            Sxx(:,kk) = mfactor(FFp,-tmp);
            kk = kk + 1;
        end
        % slope aa
        if (par.opt_slope == on & par.opt_aa == on)
            tmp = -[d0(Gx(:,par.pindx.lslope))*dSi2Cdaa; ... %.*SIL; ...
                    -d0(G(:,par.pindx.lslope))*dSi2Cdaa]; %.*SIL];
            
            Sxx(:,kk) = mfactor(FD, tmp);
        end
        % slope bb
        if (par.opt_slope == on & par.opt_bb == on)
            tmp = -bb*[d0(Gx(par.pindx.lbb))*dSi2Cdbb; ... %.*SIL; ...
                       -d0(Gx(par.pindx.lbb))*dSi2Cdbb]; %.*SIL];

            Sxx(:,kk) = mfactor(FD, tmp);
        end
        % interpp at
        if (par.opt_interpp == on & par.opt_at == on)
            tmp = at*[-dkdat.*DSI(:,par.pindx.linterpp); ...
                      (d0(dkdat) + dPFDdat)*DSI(:,par.pindx.linterpp)];
            Sxx(:,kk) = mfactor(FFp,-tmp);
            kk = kk + 1;
        end
        % interpp bt
        if (par.opt_interpp == on & par.opt_bt == on)
            tmp = -bt*[-dkdbt.*DSI(:,par.pindx.linterpp); ...
                       (d0(dkdbt) + dPFDdbt)*DSI(:,par.pindx.linterpp)];
            Sxx(:,kk) = mfactor(FFp,-tmp);
            kk = kk + 1;
        end
        % interpp aa
        if (par.opt_interpp == on & par.opt_aa == on)
            tmp = -[d0(Gx(:,par.pindx.linterpp))*dSi2Cdaa; ... %.*SIL; ...
                    -d0(G(:,par.pindx.linterpp))*dSi2Cdaa]; %.*SIL];
            
            Sxx(:,kk) = mfactor(FD, tmp);
        end
        % interpp bb
        if (par.opt_interpp == on & par.opt_bb == on)
            tmp = -bb*[d0(Gx(par.pindx.lbb))*dSi2Cdbb; ... %.*SIL; ...
                       -d0(Gx(par.pindx.lbb))*dSi2Cdbb]; %.*SIL];

            Sxx(:,kk) = mfactor(FD, tmp);
        end

        % alpha at
        if (par.opt_alpha == on & par.opt_at == on)
            tmp = at*[-dkdat.*DSI(:,par.pindx.lalpha); ...
                      (d0(dkdat) + dPFDdat)*DSI(:,par.pindx.lalpha)];
            Sxx(:,kk) = mfactor(FFp,-tmp);
            kk = kk + 1;
        end
        % alpha bt
        if (par.opt_alpha == on & par.opt_bt == on)
            tmp = -bt*[-dkdbt.*DSI(:,par.pindx.lalpha); ...
                       (d0(dkdbt) + dPFDdbt)*DSI(:,par.pindx.lalpha)];
            Sxx(:,kk) = mfactor(FFp,-tmp);
            kk = kk + 1;
        end
        % alpha aa
        if (par.opt_alpha == on & par.opt_aa == on)
            tmp = -Gxx(:,kk);
            Sxx(:,kk) = mfactor(FD, tmp);
        end
        % alpha bb
        if (par.opt_alpha == on & par.opt_bb == on)
            tmp = -Gxx(:,kk);
            Sxx(:,kk) = mfactor(FD, tmp);
        end

        % beta at
        if (par.opt_beta == on & par.opt_at == on)
            tmp = at*[-dkdat.*DSI(:,par.pindx.lbeta); ...
                      (d0(dkdat) + dPFDdat)*DSI(:,par.pindx.lbeta)];
            Sxx(:,kk) = mfactor(FFp,-tmp);
            kk = kk + 1;
        end
        % beta bt
        if (par.opt_beta == on & par.opt_bt == on)
            tmp = -bt*[-dkdbt.*DSI(:,par.pindx.lbeta); ...
                       (d0(dkdbt) + dPFDdbt)*DSI(:,par.pindx.lbeta)];
            Sxx(:,kk) = mfactor(FFp,-tmp);
            kk = kk + 1;
        end
        % beta aa
        if (par.opt_beta == on & par.opt_aa == on)
            tmp = -Gxx(:,kk);        
            Sxx(:,kk) = mfactor(FD, tmp);
        end
        % beta bb
        if (par.opt_beta == on & par.opt_bb == on)
            tmp = -Gxx(:,kk);
            Sxx(:,kk) = mfactor(FD, tmp);
        end
        % at at
        if (par.opt_at == on)
            vout = buildPFD(par, parm);
            d2PFDdat2 = vout.d2PFDdat2;
            tmp = -2*at*[-dkdat.*DSI(:,par.pindx.lat); ...
                         (d0(dkdat)+dPFDdat)*DSI(:,par.pindx.lat)] + ...
                  -at*[-dkdat.*DSI; ...
                       dkdat.*DSI + d2PFDdat2*DSI];
            
            Sxx(:,kk) = mfactor(FD, tmp);
        end
        % at bt
        if (par.opt_at == on & par.opt_bt == on)
            vout = buildPFD(par, parm);
            d2PFDdatdbt = vout.d2PFDdatdbt;
            tmp = -at*[-dkdat.*DSI(:,par.pindx.lbt); ...
                       (d0(dkdat)+dPFDdat)*DSI(:,par.pindx.lbt)] + ...
                  -bt*[-dkdbt.*DSI(:,par.pindx.lat); ...
                       (d0(dkdbt)+dPFDdbt)*DSI(:,par.pindx.lat)] + ...
                  -at*[-dkdat.*DSI; ...
                       dkdat.*DSI + d2PFDdatdbt*DSI];
            
            Sxx(:,kk) = mfactor(FD, tmp);
        end
        % at aa
        if (par.opt_at == on & par.opt_aa == on)
            tmp = -at*[-dkdat.*DSI(:,par.pindx.laa); ...
                       (d0(dkdat) + dPFDdat)*DSI(:,par.pindx.laa)] + ...
                  -[d0(Gx(:,par.pindx.lat))*dSi2Cdaa; ... %.*SIL; ...
                    -d0(Gx(:,par.pindx.lat))*dSi2Cdaa]; %.*SIL];
            
            Sxx(:,kk) = mfactor(FD, tmp);
        end
        
        % at bb
        if (par.opt_at == on & par.opt_bb == on)
            tmp = -at*[-dkdat.*DSI(:,par.pindx.lbb); ...
                       (d0(dkdat) + dPFDdat)*DSI(:,par.pindx.lbb)] + ...
                  -[d0(Gx(:,par.pindx.lat))*dSi2Cdbb; ... %.*SIL; ...
                    -d0(Gx(:,par.pindx.lat))*dSi2Cdbb]; %.*SIL];

            Sxx(:,kk) = mfactor(FD, tmp);
        end
        % bt bt
        if (par.opt_bt == on)
            vout = buildPFD(par, parm);
            d2PFDdbt2 = vout.d2PFDdbt2;
            tmp = -2*bt*[-dkdbt.*DSI(:,par.pindx.lbt); ...
                         (d0(dkdbt)+d2PFDdbt2)*DSI(:,par.pindx.lbt)] + ...
                  -bt*[-dkdbt.*DSI; ...
                       dkdbt.*DSI + d2PFDdbt2*DSI];
            
            Sxx(:,kk) = mfactor(FD, tmp);
        end
        % bt aa
        if (par.opt_bt == on & par.opt_aa == on)
            tmp = -bt*[-dkdbt.*DSI(:,par.pindx.laa); ...
                       (d0(dkdbt) + dPFDdbt)*DSI(:,par.pindx.laa)] + ...
                  -[d0(Gx(:,par.pindx.lbt))*dSi2Cdaa; ... %.*SIL; ...
                    -d0(Gx(:,par.pindx.lbt))*dSi2Cdaa]; %.*SIL];
            
            Sxx(:,kk) = mfactor(FD, tmp);
        end
        % bt bb
        if (par.opt_bt == on & par.opt_bb == on)
            tmp = -bt*[-dkdbt.*DSI(:,par.pindx.lbb); ...
                       (d0(dkdbt) + dPFDdbt)*DSI(:,par.pindx.lbb)] + ...
                  -[d0(Gx(:,par.pindx.lbt))*dSi2Cdbb; ... %.*SIL; ...
                    -d0(Gx(:,par.pindx.lbt))*dSi2Cdbb]; %.*SIL];

            Sxx(:,kk) = mfactor(FD, tmp);
        end
        % aa aa
        if (par.opt_aa == on & par.opt_aa == on)
            tmp = -Gxx(:,kk);
            Sxx(:,kk) = mfactor(FD, tmp);
        end
        % aa bb
        if (par.opt_aa == on & par.opt_bb == on)
            tmp = -Gxx(:,kk);
            Sxx(:,kk) = sparse(nwet,1);
        end
        % bb bb
        if (par.opt_bb == on & par.opt_bb == on)
            tmp = -Gxx(:,kk);
            Sxx(:,kk) = mfactor(FD, tmp);
        end    
    end 
    
end
%% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++