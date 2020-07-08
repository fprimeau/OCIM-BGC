function [f, fx, fxx] = neglogpost(x, parm, par)
global iter
on = true; off = false;
if iter <= 10
    % reset parameters if they are too large/small;
    [parm,x] = reset_par(x, parm, par);
end
fprintf('current iteration is %d \n',iter);
iter = iter + 1;
% print out current parameter values to the log file
PrintPara(x, par);
%
nx = length(x); % number of parameters
dVt  = parm.dVt  ;
M3d  = parm.M3d  ;
iwet = parm.iwet ;
nwet = parm.nwet ;
%
f = 0;
%%%%%%%%%%%%%%%%%%   Solve P    %%%%%%%%%%%%%%%%%%%%%%%%
%
idip = find(parm.po4raw(iwet)>0);
Wp   = d0(dVt(iwet(idip))/sum(dVt(iwet(idip))));
mu_dip  = sum(Wp*parm.po4raw(iwet(idip)))/sum(diag(Wp));
var_dip = sum(Wp*(parm.po4raw(iwet(idip))-mu_dip).^2)/sum(diag(Wp));
Wp = Wp/var_dip;
%
[parm, P, Px, Pxx] = eqPcycle(par, parm, x);
DIP = M3d+nan;  DIP(iwet) = P(1+0*nwet:1*nwet) ;
POP = M3d+nan;  POP(iwet) = P(1+1*nwet:2*nwet) ;
DOP = M3d+nan;  DOP(iwet) = P(1+2*nwet:3*nwet) ;
parm.Px  = Px;
parm.Pxx = Pxx;
parm.DIP = DIP(iwet);
% DIP error
ep = DIP(iwet(idip)) - parm.po4raw(iwet(idip));
f = f + 0.5*(ep.'*Wp*ep);
fname = strcat(parm.VER,'_P');
save(fname,'DIP','DOP','POP')
%%%%%%%%%%%%%%%%%%   End Solve P    %%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%   Solve Si       %%%%%%%%%%%%%%%%%%%%
if (par.Simodel == on)
    %
    isil = find(parm.sio4raw(iwet)>0);
    Ws   = d0(dVt(iwet(isil))/sum(dVt(iwet(isil))));
    mu_sil = sum(Ws*parm.sio4raw(iwet(isil)))/sum(diag(Ws));
    var_sil = sum(Ws*(parm.sio4raw(iwet(isil))-mu_sil).^2)/sum(diag(Ws));
    Ws = Ws/var_sil;
    %
    [parm,Si,Six,Sixx] = eqSicycle(par, parm, x);
    SIL = M3d+nan;  SIL(iwet) = Si(1:nwet);
    DSI = M3d+nan;  DSI(iwet) = Si(nwet+1:end);
    % SiO error
    es = SIL(iwet(isil)) - parm.sio4raw(iwet(isil));
    f = f + 0.5*(es.'*Ws*es);
    fname = strcat(parm.VER,'_Si');
    save(fname,'SIL','DSI')
end
%%%%%%%%%%%%%%%%%%   End Solve Si    %%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%     Solve C   %%%%%%%%%%%%%%%%%%%%%%%%
if (par.Cmodel == on)
    idic = find(parm.dicraw(iwet)>0);
    Wc   = d0(dVt(iwet(idic))/sum(dVt(iwet(idic))));
    mu_dic = sum(Wc*parm.dicraw(iwet(idic)))/sum(diag(Wc));
    var_dic = sum(Wc*(parm.dicraw(iwet(idic))-mu_dic).^2)/sum(diag(Wc));
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
    ec  = DIC(iwet(idic)) - parm.dicraw(iwet(idic));
    f   = f + 0.5*(ec.'*Wc*ec);
    %
    fname = strcat(parm.VER,'_C');
    save(fname,'DIC', 'POC', 'DOC', 'CaC')
end
%%%%%%%%%%%%%%%%%%   End Solve C    %%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%   Solve O    %%%%%%%%%%%%%%%%%%%%%%%%
if (par.Omodel == on)
    %
    io2 = find(parm.o2raw(iwet)>0);
    Wo  = d0(dVt(iwet(io2))/sum(dVt(iwet(io2))));
    mu_o2 = sum(Wo*parm.o2raw(iwet(io2)))/sum(diag(Wo));
    var_o2 = sum(Wo*(parm.o2raw(iwet(io2))-mu_o2).^2)/sum(diag(Wo));
    Wo = Wo/var_o2;
    %
    [parm, O, Ox, Oxx] = eqOcycle(par, parm, x);
    O2 = M3d+nan;  O2(iwet) = O;
    %
    eo  = O2(iwet(io2)) - parm.o2raw(iwet(io2));
    f   = f + 0.5*(eo.'*Wo*eo);
    %
    fname = strcat(parm.VER,'_O2');
    save(fname,'O2')
end
%%%%%%%%%%%%%%%%%%   End Solve O    %%%%%%%%%%%%%%%%%%%%
fprintf('current objective function value is %3.3e \n',f);

%% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
if (nargout > 1)
    fx = zeros(length(x),1);
    px  = Px(1:nwet,:);
    npx = parm.npx;
    % ---------------------------------
    for ji = 1 : npx
        fx(ji) = ep.'*Wp*px(idip,ji);
    end
    % ---------------------------------
    if (par.Simodel == on)
        sx  = Six(1:nwet,:);
        nsx = parm.nsx;
        % --------------------------
        for ji = 1 : npx+nsx
            fx(ji) = fx(ji) + es.'*Ws*sx(isil, ji);
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
            fx(ji) = fx(ji) + eo.'*Wo*ox(io2, ji);
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
            fxx(ju,jo) = px(idip,ju).'*Wp*px(idip,jo) + ep.'*Wp*pxx(idip,kk);
            % Simodel
            if par.Simodel == on
                sxx = Sixx(1:nwet,:);            
                fxx(ju,jo) = fxx(ju, jo) + ...
                    sx(isil,ju).'*Ws*sx(isil,jo) + es.'*Ws*sxx(isil,kk);
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
                    ox(io2,ju).'*Wo*ox(io2,jo) + eo.'*Wo*oxx(io2,kk);
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
                    sx(isil,ju).'*Ws*sx(isil,jo) +es.'*Ws*sxx(isil,kk);
                
                fxx(jo, ju) = fxx(ju, jo);
            end 
        end
        % -----------------------------
        % Only Si parameters
        for ju = (npx+1):(npx+nsx)
            for jo = ju:(npx+nsx)
                kk = kk + 1;
                fxx(ju, jo) = fxx(ju, jo) + ...
                    sx(isil,ju).'*Ws*sx(isil,jo) + es.'*Ws*sxx(isil,kk);
                
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
                    ox(io2,ju).'*Wo*ox(io2,jo) + eo.'*Wo*oxx(io2,kk);
                
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
                    ox(io2,ju).'*Wo*ox(io2,jo) + eo.'*Wo*oxx(io2,kk);
                
                fxx(jo, ju) = fxx(ju, jo);
            end 
        end
        % -------------------------------
        % P and O model parameters
        for ju = 1:npx
            for jo = (npx+ncx+1):(npx+ncx+nox)
                kk = kk + 1;
                fxx(ju, jo) = fxx(ju, jo) + ...
                    ox(io2,ju).'*Wo*ox(io2,jo) + eo.'*Wo*oxx(io2,kk);
                
                fxx(jo, ju) = fxx(ju, jo);
            end 
        end
        % -------------------------------
        % C and O model parameters
        for ju = (npx+1):(npx+ncx)
            for jo = (npx+ncx+1):(npx+ncx+nox)
                kk = kk + 1;
                fxx(ju, jo) = fxx(ju, jo) + ...
                    ox(io2,ju).'*Wo*ox(io2,jo) + eo.'*Wo*oxx(io2,kk);
                
                fxx(jo, ju) = fxx(ju, jo);
            end 
        end
        % -------------------------------
        % O model parameters
        for ju = (npx+ncx+1):(npx+ncx+nox)
            for jo = ju:(npx+ncx+nox)
                kk = kk + 1;
                fxx(ju, jo) = fxx(ju, jo) + ...
                    ox(io2,ju).'*Wo*ox(io2,jo) + eo.'*Wo*oxx(io2,kk);
                
                fxx(jo, ju) = fxx(ju, jo);
            end 
        end
    end
end

