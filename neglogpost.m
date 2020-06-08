function [f, fx, fxx] = neglogpost(x, parm, par)
on = true; off = false;
% reset parameters if they are too large/small;
% [x] = reset_par(x, parm, par);
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
parm.Pxx = Pxx;
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
    [parm,Si,Six,Sixx] = eqSicycle(par, parm, x);
    SIO = M3d+nan;  SIO(iwet) = Si(1:nwet);
    DSI = M3d+nan;  DSI(iwet) = Si(nwet+1:end);
    % SiO error
    es = SIO(iwet) - parm.SIL(iwet);
    f = f + 0.5*(es.'*Ws*es); 
end
%%%%%%%%%%%%%%%%%%   End Solve Si    %%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%%%%%%%%%%%     Slove C   %%%%%%%%%%%%%%%%%%%%%%%%
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
    %
    save tmpC DIC POC DOC CaC
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
    %
    eo  = O - parm.o2obs;
    f   = f + 0.5*(eo.'*Wo*eo);
    %
    save tmpO O2
end
%%%%%%%%%%%%%%%%%%   End Solve O    %%%%%%%%%%%%%%%%%%%%
fprintf('current objective function value is %3.3e \n',f);
%
%% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
if (nargout > 1)
    fx = zeros(length(x),1);
    parm_index = [];
    if par.opt_sigma == on 
        parm_index = [parm_index, par.pindx.lsigma];
    end
    if par.opt_kappa_dp == on
        parm_index = [parm_index, par.pindx.lkappa_dp];
    end
    if par.opt_slopep == on
        parm_index = [parm_index, par.pindx.slopep];
    end 
    if par.opt_interpp == on 
        parm_index = [parm_index, par.pindx.linterpp];
    end
    if par.opt_alpha == on
        parm_index = [parm_index, par.pindx.lalpha];
    end
    if par.opt_beta == on
        parm_index = [parm_index, par.pindx.lbeta];
    end
    % ---------------------------------
    for ji = 1:length(parm_index)
        fx(ji) = ep.'*Wp*Px(1:nwet,ji);
    end
    % ---------------------------------
    %
    if (nargout>2)
        npx = parm.npx;    
        fxx = [];
        px  = Px(1:nwet,:);
        pxx = Pxx(1:nwet,:);
        % 
        kk = 1;
        % ----------------------------------------------------------------
        for ju = 1:npx
            for jo = ju:npx
                fxx(parm_index(ju),parm_index(jo)) = ...
                    px(:,parm_index(ju)).'*Wp*px(:,parm_index(jo)) + ...
                    ep.'*Wp*pxx(:,kk);
                if par.Simodel == on
                    sx  = Six(1:nwet,:);
                    sxx = Sixx(1:nwet,:);            
                    fxx(parm_index(ju),parm_index(jo)) = ...
                        fxx(parm_index(ju), parm_index(jo)) + ...
                        sx(:,parm_index(ju)).'*Ws*sx(:,parm_index(jo)) + ...
                        es.'*Ws*sxx(:,kk);
                end 
                fxx(parm_index(jo),parm_index(ju)) = ...
                    fxx(parm_index(ju),parm_index(jo));
                kk = kk + 1;
            end 
        end
        % ----------------------------------------------------------------
    end
    
    if (par.Simodel == on)
        nsx = parm.nsx;
        if par.opt_bsi == on
            parm_index = [parm_index, par.pindx.lbsi];
        end 
        if par.opt_at == on
            parm_index = [parm_index, par.pindx.lat];
        end 
        if par.opt_bt == on
            parm_index = [parm_index, par.pindx.lbt];
        end 
        if par.opt_aa == on
            parm_index = [parm_index, par.pindx.aa];
        end 
        if par.opt_bb == on 
            parm_index = [parm_index, par.pindx.lbb];
        end
        if par.opt_kappa_gs == on
            parm_index = [parm_index, par.pindx.lkappa_gs];
        end
        % --------------------------
        for ji = 1:length(parm_index)
            fx(ji) = fx(ji) + es.'*Ws*sx(:, ji);
        end
        % --------------------------
        %
        if (nargout >2)
            tmp = sparse(npx+nsx, npx+nsx);
            tmp(1:npx,1:npx) = fxx;
            fxx = tmp;
            % ------------------------------------------------------------
            for ju = 1:npx
                for jo = (npx+1):(npx+nsx)
                    fxx(parm_index(ju), parm_index(jo)) = ...
                        fxx(parm_index(ju), parm_index(jo)) + ...
                        sx(:,parm_index(ju)).'*Ws*sx(:,parm_index(jo)) + ...
                        es.'*Ws*sxx(:,kk);
                    
                    fxx(parm_index(jo), parm_index(ju)) = ...
                        fxx(parm_index(ju), parm_index(jo));
                    kk = kk + 1;
                end 
            end
            % ------------------------------------------------------------
            for ju = (npx+1):(npx+nsx)
                for jo = ju:(npx+nsx)
                    fxx(parm_index(ju), parm_index(jo)) = ...
                        fxx(parm_index(ju), parm_index(jo)) + ...
                        sx(:,parm_index(ju)).'*Ws*sx(:,parm_index(jo)) + ...
                        es.'*Ws*sxx(:,kk);
                    
                    fxx(parm_index(jo), parm_index(ju)) = ...
                        fxx(parm_index(ju), parm_index(jo));
                    kk = kk + 1;
                end 
            end
        end
        % ----------------------------------------------------------------
    end
end 
