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
EXP = 'SOxhalf';
%
f = 0;
%%%%%%%%%%%%%%%%%%   Solve P    %%%%%%%%%%%%%%%%%%%%%%%%
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
fname = strcat(EXP,'_P');
save(fname,'DIP','DOP','POP')
%%%%%%%%%%%%%%%%%%   End Solve P    %%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%%%%%%%%%%%   Solve Si       %%%%%%%%%%%%%%%%%%%%
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
%%%%%%%%%%%%%%%%%%     Solve C   %%%%%%%%%%%%%%%%%%%%%%%%
if (par.Cmodel == on)
    mu_dic = sum(W*parm.DICobs(iwet))/sum(diag(W));
    var_dic = sum(W*(parm.DICobs(iwet)-mu_dic).^2)/sum(diag(W));
    Wc = W/var_dic;
    
    [parm, C, Cx, Cxx] = eqCcycle(par, parm, x);
    DIC = M3d+nan; DIC(iwet) = C(0*nwet+1:1*nwet) ;
    POC = M3d+nan; POC(iwet) = C(1*nwet+1:2*nwet) ;
    DOC = M3d+nan; DOC(iwet) = C(2*nwet+1:3*nwet) ;
    CaC = M3d+nan; CaC(iwet) = C(3*nwet+1:4*nwet) ;
    parm.DIC = DIC(iwet);      parm.DOC = DOC(iwet);
    parm.DICx = Cx(1:nwet,:);  parm.DOCx = Cx(2*nwet+1:3*nwet,:);
    parm.DICxx = Cxx(1:nwet,:);  parm.DOCxx = Cxx(2*nwet+1:3*nwet,:);
    % DIC error
    DIC = DIC + parm.human_co2;
    ec  = DIC(iwet) - parm.DICobs(iwet);
    f   = f + 0.5*(ec.'*Wc*ec);
    %
    fname = strcat(EXP,'_C');
    save(fname,'DIC', 'POC', 'DOC', 'CaC')
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
    [parm, O, Ox, Oxx] = eqOcycle(par, parm, x);
    O2 = M3d+nan;  O2(iwet) = O;
    %
    eo  = O - parm.o2obs;
    f   = f + 0.5*(eo.'*Wo*eo);
    %
    fname = strcat(EXP,'_O2');
    save(fname,'O2')
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
                % Simodel
                if par.Simodel == on
                    sx  = Six(1:nwet,:);
                    sxx = Sixx(1:nwet,:);            
                    fxx(parm_index(ju),parm_index(jo)) = ...
                        fxx(parm_index(ju), parm_index(jo)) + ...
                        sx(:,parm_index(ju)).'*Ws*sx(:,parm_index(jo)) + ...
                        es.'*Ws*sxx(:,kk);
                end 
                % Cmodel
                if par.Cmodel == on
                    cx  = Cx(1:nwet,:);
                    cxx = Cxx(1:nwet,:);            
                    fxx(parm_index(ju),parm_index(jo)) = ...
                        fxx(parm_index(ju), parm_index(jo)) + ...
                        cx(:,parm_index(ju)).'*Wc*cx(:,parm_index(jo)) + ...
                        ec.'*Wc*cxx(:,kk);
                end
                % Omodel
                if par.Omodel == on
                    ox  = Ox(1:nwet,:);
                    oxx = Oxx(1:nwet,:);            
                    fxx(parm_index(ju),parm_index(jo)) = ...
                        fxx(parm_index(ju), parm_index(jo)) + ...
                        ox(:,parm_index(ju)).'*Wo*ox(:,parm_index(jo)) + ...
                        eo.'*Wo*oxx(:,kk);
                end 
                % make Hessian symetric;
                fxx(parm_index(jo),parm_index(ju)) = ...
                    fxx(parm_index(ju),parm_index(jo));
                kk = kk + 1;
            end 
        end
        % ----------------------------------------------------------------
    end
    
    if (par.Simodel == on)
        sx  = Six(1:nwet,:);
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
            % P model parameters and Si parameters
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
            % Only Si parameters
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
    % C model
    if (par.Cmodel == on & par.Omodel == off)
        cx  = Cx(1:nwet,:);
        ncx = parm.ncx;
        % --------------------------
        for ji = 1:(npx+ncx)
            fx(ji) = fx(ji) + ec.'*Wc*cx(:, ji);
        end
        % --------------------------
        %
        if (nargout >2)
            tmp = sparse(npx+ncx, npx+ncx);
            tmp(1:npx,1:npx) = fxx;
            fxx = tmp;
            % ------------------------------------------------------------
            % P model parameters and C parameters
            for ju = 1:npx
                for jo = (npx+1):(npx+ncx)
                    fxx(ju, jo) = ...
                        fxx(ju, jo) + ...
                        cx(:,ju).'*Wc*cx(:,jo) + ec.'*Wc*cxx(:,kk);
                    fxx(jo,ju) = fxx(ju, jo);
                    kk = kk + 1;
                end 
            end
            % ------------------------------------------------------------
            % Only C parameters
            for ju = (npx+1):(npx+ncx)
                for jo = ju:(npx+ncx)
                    fxx(ju, jo) = fxx(ju, jo) + ...
                        cx(:,ju).'*Wc*cx(:,jo) + ec.'*Wc*cxx(:,kk);
                    fxx(jo, ju) = fxx(ju, jo);
                    kk = kk + 1;
                end 
            end
        end
        % ----------------------------------------------------------------
    end
    % C model and O model are both on
    if (par.Cmodel == on & par.Omodel == on)
        cx  = Cx(1:nwet,:);
        ox  = Ox(1:nwet,:);
        nox = parm.nox;
        ncx = parm.ncx;
        % --------------------------
        % add C model Gradient
        for ji = 1:(npx+ncx)
            fx(ji) = fx(ji) + ec.'*Wc*cx(:, ji);
        end
        % add O model Gradient
        for ji = 1:(npx+ncx+nox)
            fx(ji) = fx(ji) + eo.'*Wo*ox(:, ji);
        end
        % --------------------------
        %
        if (nargout > 2)
            tmp = sparse(npx+ncx+nox, npx+ncx+nox);
            tmp(1:npx,1:npx) = fxx;
            fxx = tmp;
            % ------------------------------------------------------------
            % P and C model parameters
            for ju = 1:npx
                for jo = (npx+1):(npx+ncx)
                    fxx(ju, jo) = ...
                        fxx(ju, jo) + ...
                        cx(:,ju).'*Wc*cx(:,jo) + ec.'*Wc*cxx(:,kk) + ...
                        ox(:,ju).'*Wo*ox(:,jo) + eo.'*Wo*oxx(:,kk);
                    
                    fxx(jo, ju) = fxx(ju, jo);
                    kk = kk + 1;
                end 
            end
            % ------------------------------------------------------------
            % C model parameters
            for ju = (npx+1):(npx+ncx)
                for jo = ju:(npx+ncx)
                    fxx(ju, jo) = ...
                        fxx(ju, jo) + ...
                        cx(:,ju).'*Wc*cx(:,jo) + ec.'*Wc*cxx(:,kk) + ...
                        ox(:,ju).'*Wo*ox(:,jo) + eo.'*Wo*oxx(:,kk);
   
                    fxx(jo, ju) = fxx(ju, jo);
                    kk = kk + 1;
                end 
            end
            % ------------------------------------------------------------
            % P and O model parameters
            for ju = 1:npx
                for jo = (npx+ncx+1):(npx+ncx+nox)
                    fxx(ju, jo) = ...
                        fxx(ju, jo) + ...
                        ox(:,ju).'*Wo*ox(:,jo) + eo.'*Wo*oxx(:,kk);

                    fxx(jo, ju) = fxx(ju, jo);
                    kk = kk + 1;
                end 
            end
            % ------------------------------------------------------------
            % C and O model parameters
            for ju = (npx+1):(npx+ncx)
                for jo = (npx+ncx+1):(npx+ncx+nox)
                    fxx(ju, jo) = ...
                        fxx(ju, jo) + ...
                        ox(:,ju).'*Wo*ox(:,jo) + eo.'*Wo*oxx(:,kk);
                    
                    fxx(jo, ju) = fxx(ju, jo);
                    kk = kk + 1;
                end 
            end
            % ------------------------------------------------------------
            % O model parameters
            for ju = (npx+ncx+1):(npx+ncx+nox)
                for jo = ju:(npx+ncx+nox)
                    fxx(ju, jo) = ...
                        fxx(ju, jo) + ...
                        ox(:,ju).'*Wo*ox(:,jo) + eo.'*Wo*oxx(:,kk);

                    fxx(jo, ju) = fxx(ju, jo);
                    kk = kk + 1;
                end 
            end
        end
        % ----------------------------------------------------------------
    end
end 
