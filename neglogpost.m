function [f, fx, fxx] = neglogpost(x, parm, par)
on = true; off = false;
% reset parameters if they are too large/small;
[parm,x] = reset_par(x, parm, par);
nx = length(x); % number of parameters

dVt  = parm.dVt  ;
M3d  = parm.M3d  ;
iwet = parm.iwet ;
nwet = parm.nwet ;
%
EXP = 'temp_dep_b';
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
% save(fname,'DIP','DOP','POP')
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
    idic = find(parm.DICobs(iwet)>0);
    Wc   = d0(dVt(iwet(idic))/sum(dVt(iwet(idic))));
    mu_dic = sum(Wc*parm.DICobs(iwet(idic)))/sum(diag(Wc));
    var_dic = sum(Wc*(parm.DICobs(iwet(idic))-mu_dic).^2)/sum(diag(Wc));
    Wc = Wc/var_dic;
    
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
    ec  = DIC(iwet(idic)) - parm.DICobs(iwet(idic));
    f   = f + 0.5*(ec.'*Wc*ec);
    %
    fname = strcat(EXP,'_C');
    % save(fname,'DIC', 'POC', 'DOC', 'CaC')
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
    % save(fname,'O2')
end
%%%%%%%%%%%%%%%%%%   End Solve O    %%%%%%%%%%%%%%%%%%%%
fprintf('current objective function value is %3.3e \n',f);
%
%% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
if (nargout > 1)
    fx = zeros(length(x),1);
    px  = Px(1:nwet,:);
    npx = parm.npx;
    % ---------------------------------
    for ji = 1 : npx
        fx(ji) = ep.'*Wp*px(:,ji);
    end
    % ---------------------------------
    if (par.Simodel == on)
        sx  = Six(1:nwet,:);
        nsx = parm.nsx;
        % --------------------------
        for ji = 1 : npx+nsx
            fx(ji) = fx(ji) + es.'*Ws*sx(:, ji);
        end
    end 
    % ----------------------------------
    % ----------------------------------
    if (par.Cmodel == on)
        cx  = Cx(1:nwet,:);
        ncx = parm.ncx;
        % --------------------------
        for ji = 1 : npx+ncx
            fx(ji) = fx(ji) + ec.'*Wc*cx(idic, ji);
        end
    end 
    % ----------------------------------
    % ----------------------------------
    if (par.Omodel == on)
        ox  = Ox(1:nwet,:);
        nox = parm.nox;
        % --------------------------
        for ji = 1 : npx+ncx+nox
            fx(ji) = fx(ji) + eo.'*Wo*ox(:, ji);
        end
    end 
    % ----------------------------------
end 

%
if (nargout>2)
    fxx = sparse(npx,npx);
    pxx = Pxx(1:nwet,:);
    % ----------------------------------------------------------------
    kk = 0;
    for ju = 1:npx
        for jo = ju:npx
            kk = kk + 1;
            fxx(ju,jo) = px(:,ju).'*Wp*px(:,jo) + ep.'*Wp*pxx(:,kk);
            % Simodel
            if par.Simodel == on
                sxx = Sixx(1:nwet,:);            
                fxx(ju,jo) = fxx(ju, jo) + ...
                    sx(:,ju).'*Ws*sx(:,jo) + es.'*Ws*sxx(:,kk);
            end 
            % Cmodel
            if par.Cmodel == on
                cxx = Cxx(1:nwet,:);            
                fxx(ju,jo) = fxx(ju, jo) + ...
                    cx(idic,ju).'*Wc*cx(idic,jo) + ec.'*Wc*cxx(idic,kk);
            end
            % Omodel
            if par.Omodel == on
                oxx = Oxx(1:nwet,:);            
                fxx(ju,jo) = fxx(ju, jo) + ...
                    ox(:,ju).'*Wo*ox(:,jo) + eo.'*Wo*oxx(:,kk);
            end 
            % make Hessian symetric;
            fxx(jo, ju) = fxx(ju, jo);  
        end 
    end
    % ----------------------------------------------------------------
    % Si model
    if (par.Simodel == on)
        % --------------------------
        tmp = sparse(npx+nsx, npx+nsx);
        tmp(1:npx,1:npx) = fxx;
        fxx = tmp;
        % --------------------------
        % P model parameters and Si parameters
        for ju = 1:npx
            for jo = (npx+1):(npx+nsx)
                kk = kk + 1;
                fxx(ju, jo) = fxx(ju, jo) + ...
                    sx(:,ju).'*Ws*sx(:,jo) +es.'*Ws*sxx(:,kk);
                
                fxx(jo, ju) = fxx(ju, jo);
            end 
        end
        % -----------------------------
        % Only Si parameters
        for ju = (npx+1):(npx+nsx)
            for jo = ju:(npx+nsx)
                kk = kk + 1;
                fxx(ju, jo) = fxx(ju, jo) + ...
                    sx(:,ju).'*Ws*sx(:,jo) + es.'*Ws*sxx(:,kk);
                
                fxx(jo, ju) = fxx(ju, jo);
            end 
        end  
    end
    % ----------------------------------------------------------------
    % C model
    if (par.Cmodel == on & par.Omodel == off)
        %
        tmp = sparse(npx+ncx, npx+ncx);
        tmp(1:npx,1:npx) = fxx;
        fxx = tmp;
        % -------------------------------
        % P model parameters and C parameters
        for ju = 1:npx
            for jo = (npx+1):(npx+ncx)
                kk = kk + 1;
                fxx(ju, jo) = fxx(ju, jo) + ...
                    cx(idic,ju).'*Wc*cx(idic,jo) + ec.'*Wc*cxx(idic,kk);
                fxx(jo,ju) = fxx(ju, jo);
            end 
        end
        % -------------------------------
        % Only C parameters
        for ju = (npx+1):(npx+ncx)
            for jo = ju:(npx+ncx)
                kk = kk + 1;
                fxx(ju, jo) = fxx(ju, jo) + ...
                    cx(idic,ju).'*Wc*cx(idic,jo) + ec.'*Wc*cxx(idic,kk);
                fxx(jo, ju) = fxx(ju, jo);
            end 
        end
    end
    % ----------------------------------------------------------------
    % C and O model
    if (par.Cmodel == on & par.Omodel == on)
        %
        tmp = sparse(npx+ncx+nox, npx+ncx+nox);
        tmp(1:npx,1:npx) = fxx;
        fxx = tmp;
        % -------------------------------
        % P and C model parameters
        for ju = 1:npx
            for jo = (npx+1):(npx+ncx)
                kk = kk + 1;
                fxx(ju, jo) = ...
                    fxx(ju, jo) + ...
                    cx(idic,ju).'*Wc*cx(idic,jo) + ec.'*Wc*cxx(idic,kk) + ...
                    ox(:,ju).'*Wo*ox(:,jo) + eo.'*Wo*oxx(:,kk);
                
                fxx(jo, ju) = fxx(ju, jo);
            end 
        end
        % -------------------------------
        % C model parameters
        for ju = (npx+1):(npx+ncx)
            for jo = ju:(npx+ncx)
                kk = kk + 1;
                fxx(ju, jo) = ...
                    fxx(ju, jo) + ...
                    cx(idic,ju).'*Wc*cx(idic,jo) + ec.'*Wc*cxx(idic,kk) + ...
                    ox(:,ju).'*Wo*ox(:,jo) + eo.'*Wo*oxx(:,kk);
                
                fxx(jo, ju) = fxx(ju, jo);
            end 
        end
        % -------------------------------
        % P and O model parameters
        for ju = 1:npx
            for jo = (npx+ncx+1):(npx+ncx+nox)
                kk = kk + 1;
                fxx(ju, jo) = fxx(ju, jo) + ...
                    ox(:,ju).'*Wo*ox(:,jo) + eo.'*Wo*oxx(:,kk);
                
                fxx(jo, ju) = fxx(ju, jo);
            end 
        end
        % -------------------------------
        % C and O model parameters
        for ju = (npx+1):(npx+ncx)
            for jo = (npx+ncx+1):(npx+ncx+nox)
                kk = kk + 1;
                fxx(ju, jo) = fxx(ju, jo) + ...
                    ox(:,ju).'*Wo*ox(:,jo) + eo.'*Wo*oxx(:,kk);
                
                fxx(jo, ju) = fxx(ju, jo);
            end 
        end
        % -------------------------------
        % O model parameters
        for ju = (npx+ncx+1):(npx+ncx+nox)
            for jo = ju:(npx+ncx+nox)
                kk = kk + 1;
                fxx(ju, jo) = fxx(ju, jo) + ...
                    ox(:,ju).'*Wo*ox(:,jo) + eo.'*Wo*oxx(:,kk);
                
                fxx(jo, ju) = fxx(ju, jo);
            end 
        end
    end
end

