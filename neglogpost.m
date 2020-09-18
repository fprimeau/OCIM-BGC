function [f, fx, fxx, data] = neglogpost(x, par)
    global iter
    on = true; off = false;
    % print and save current parameter values to
    % a file that is used to reset parameters ;
    if iter == 0
        PrintPara(x, par) ;
    end
    % reset parameters if optimization routine
    % suggests strange parameter values ;
    if iter < 10
        x = ResetPara(x, par) ;
    end
    % print current parameters 
    if iter > 0
        PrintPara(x, par) ;    
    end
    fprintf('current iteration is %d \n',iter) ;
    iter = iter + 1  ;
    
    nx   = length(x) ; % number of parameters
    dVt  = par.dVt   ;
    M3d  = par.M3d   ;
    iwet = par.iwet  ;
    nwet = par.nwet  ;
    %
    f    = 0 ;
    %%%%%%%%%%%%%%%%%%   Solve P    %%%%%%%%%%%%%%%%%%%%%%%%
    idip = find(par.po4raw(iwet)>0) ;
    Wp   = d0(dVt(iwet(idip))/sum(dVt(iwet(idip)))) ;
    mu   = sum(Wp*par.po4raw(iwet(idip)))/sum(diag(Wp)) ;
    var  = sum(Wp*(par.po4raw(iwet(idip))-mu).^2)/sum(diag(Wp)) ;
    Wp   = Wp/var ;
    %
    [par, P, Px, Pxx] = eqPcycle(x, par) ;
    DIP = M3d+nan ;  DIP(iwet) = P(1+0*nwet:1*nwet) ;
    POP = M3d+nan ;  POP(iwet) = P(1+1*nwet:2*nwet) ;
    DOP = M3d+nan ;  DOP(iwet) = P(1+2*nwet:3*nwet) ;
    par.Px   = Px  ;
    par.Pxx  = Pxx ;
    par.DIP  = DIP(iwet) ;
    data.DIP = DIP ; data.POP = POP ; data.DOP = DOP ;
    % DIP error
    ep = DIP(iwet(idip)) - par.po4raw(iwet(idip)) ;
    f  = f + 0.5*(ep.'*Wp*ep) ;
    %%%%%%%%%%%%%%%%%%   End Solve P    %%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%   Solve Si       %%%%%%%%%%%%%%%%%%%%
    if (par.Simodel == on)
        isil = find(par.sio4raw(iwet)>0) ;
        Ws   = d0(dVt(iwet(isil))/sum(dVt(iwet(isil)))) ;
        mu   = sum(Ws*par.sio4raw(iwet(isil)))/sum(diag(Ws)) ;
        var  = sum(Ws*(par.sio4raw(iwet(isil))-mu).^2)/sum(diag(Ws));
        Ws   = Ws/var ;
        %
        [par,Si,Six,Sixx] = eqSicycle(x, par)   ;
        SIL = M3d+nan ;  SIL(iwet) = Si(1:nwet) ;
        DSI = M3d+nan ;  DSI(iwet) = Si(nwet+1:end) ;
        data.SIL = SI ;  data.DSI  = DSI ;  
        % SiO error
        es = SIL(iwet(isil)) - par.sio4raw(iwet(isil)) ;
        f  = f + 0.5*(es.'*Ws*es) ;
    end
    %%%%%%%%%%%%%%%%%%   End Solve Si    %%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%     Solve C   %%%%%%%%%%%%%%%%%%%%%%%%
    if (par.Cmodel == on)
        idic = find(par.dicraw(iwet)>0) ;
        Wc   = d0(dVt(iwet(idic))/sum(dVt(iwet(idic)))) ;
        mu   = sum(Wc*par.dicraw(iwet(idic)))/sum(diag(Wc)) ;
        var  = sum(Wc*(par.dicraw(iwet(idic))-mu).^2)/sum(diag(Wc));
        Wc   = Wc/var ;
        
        [par, C, Cx, Cxx] = eqCcycle(x, par) ;
        DIC = M3d+nan ;  DIC(iwet) = C(0*nwet+1:1*nwet) ;
        POC = M3d+nan ;  POC(iwet) = C(1*nwet+1:2*nwet) ;
        DOC = M3d+nan ;  DOC(iwet) = C(2*nwet+1:3*nwet) ;
        CaC = M3d+nan ;  CaC(iwet) = C(3*nwet+1:4*nwet) ;
        par.DIC  = DIC(iwet) ;
        par.DOC  = DOC(iwet) ;
        DIC = DIC + par.human_co2  ;
        par.Cx   = Cx        ;  par.Cxx  = Cxx ;
        data.DIC = DIC       ;  data.POC = POC ;
        data.DOC = DOC       ;  data.CaC = CaC ;
        % DIC error
        ec  = DIC(iwet(idic)) - par.dicraw(iwet(idic)) ;
        f   = f + 0.5*(ec.'*Wc*ec) ;
    end
    %%%%%%%%%%%%%%%%%%   End Solve C    %%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%   Solve O    %%%%%%%%%%%%%%%%%%%%%%%%
    if (par.Omodel == on)
        io2 = find(par.o2raw(iwet)>0) ;
        Wo  = d0(dVt(iwet(io2))/sum(dVt(iwet(io2)))) ;
        mu  = sum(Wo*par.o2raw(iwet(io2)))/sum(diag(Wo)) ;
        var = sum(Wo*(par.o2raw(iwet(io2))-mu).^2)/sum(diag(Wo)) ;
        Wo  = Wo/var ;
        %
        [par, O, Ox, Oxx] = eqOcycle(x, par) ;
        O2 = M3d+nan ;  O2(iwet) = O ;
        data.O2 = O2 ;
        eo = O2(iwet(io2)) - par.o2raw(iwet(io2)) ;
        f  = f + 0.5*(eo.'*Wo*eo)   ;
    end
    %%%%%%%%%%%%%%%%%%   End Solve O    %%%%%%%%%%%%%%%%%%%%
    fprintf('current objective function value is %3.3e \n\n',f) 

    %% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    % calculate gradient
    if (nargout > 1)
        fx = zeros(length(x), 1) ;
        px  = Px(1:nwet,:)       ;
        npx = par.npx            ;
        % ---------------------------------
        for ji = 1 : npx
            fx(ji) = ep.'*Wp*px(idip,ji) ;
        end
        % ---------------------------------
        if (par.Simodel == on & par.Omodel == off & par.Cmodel == off)
            sx  = Six(1:nwet,:) ;
            nsx = par.nsx       ;
            % --------------------------
            for ji = 1 : npx+nsx
                fx(ji) = fx(ji) + es.'*Ws*sx(isil, ji) ;
            end
        end 
        % ----------------------------------
        % ----------------------------------
        if (par.Cmodel == on & par.Simodel == off)
            cx  = Cx(1:nwet, :) ;
            ncx = par.ncx       ;
            % --------------------------
            for ji = 1 : npx+ncx
                fx(ji) = fx(ji) + ec.'*Wc*cx(idic, ji) ;
            end
        end 
        % ----------------------------------
        % ----------------------------------
        if (par.Omodel == on & par.Simodel == off)
            ox  = Ox(1:nwet, :) ;
            nox = par.nox       ;
            % --------------------------
            for ji = 1 : npx+ncx+nox
                fx(ji) = fx(ji) + eo.'*Wo*ox(io2, ji) ;
            end
        end 
        % ----------------------------------
        % ----------------------------------
        if (par.Cmodel == on & par.Simodel == on & par.Omodel == off)
            ncx = par.ncx ;
            nsx = par.nsx ;
            cx  = Cx(1:nwet, :)  ;
            sx  = Six(1:nwet, :) ;
            % --------------------------
            for ji = 1 : npx
                fx(ji) = fx(ji) + ec.'*Wc*cx(idic, ji) ...
                         + es.'*Ws*sx(isil, ji) ;
            end
            % --------------------------
            for ji = npx+1 : npx+ncx
                fx(ji) = fx(ji) + ec.'*Wc*cx(idic, ji) ; 
            end
            % --------------------------
            for ji = npx+ncx+1 : npx+ncx+nsx
                fx(ji) = fx(ji) + es.'*Ws*sx(isil, ji) ;
            end
        end 
        % ----------------------------------
        % ----------------------------------
        if (par.Cmodel == on & par.Simodel == on & par.Omodel == on)
            ncx = par.ncx ;
            nox = par.nox ;
            nsx = par.nsx ;
            cx  = Cx(1:nwet,:) ;
            ox  = Ox(1:nwet,:) ;
            sx  = Six(1:nwet,:);
            % --------------------------
            for ji = 1 : npx
                fx(ji) = fx(ji) ...
                         + ec.'*Wc*cx(idic, ji) ...
                         + eo.'*Wo*ox(io2 , ji) ...
                         + es.'*Ws*sx(isil, ji) ;
            end
            % --------------------------
            for ji = npx+1 : npx+ncx
                fx(ji) = fx(ji) + ec.'*Wc*cx(idic, ji) ;
            end
            % --------------------------
            for ji = npx+1 : npx+ncx+nox
                fx(ji) = fx(ji) + eo.'*Wo*ox(io2, ji) ;
            end
            % --------------------------
            for ji = npx+ncx+nox+1 : npx+ncx+nox+nsx
                fx(ji) = fx(ji) + es.'*Ws*sx(isil, ji) ;
            end
        end 
        % ----------------------------------
    end 
    %
    if (nargout>2)
        fxx = sparse(npx, npx) ;
        pxx = Pxx(1:nwet, :)   ;
        % ----------------------------------------------------------------
        kk = 0;
        for ju = 1:npx
            for jo = ju:npx
                kk = kk + 1 ;
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
        kpp = kk;
        % ----------------------------------------------------------------
        % C model
        if (par.Cmodel == on & par.Omodel == off & par.Simodel == off)
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
        if (par.Cmodel == on & par.Omodel == on & par.Simodel == off)
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
                        ox(io2, ju).'*Wo*ox(io2, jo) + eo.'*Wo*oxx(io2, kk);
                    
                    fxx(jo, ju) = fxx(ju, jo);
                end 
            end
        end
        % ----------------------------------------------------------------
        % Si model on; Cmodel off; Omodel off;
        if (par.Cmodel == off & par.Omodel == off & par.Simodel == on)
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
                        sx(isil,ju).'*Ws*sx(isil,jo) + es.'*Ws*sxx(isil,kk);
                    
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
        % Cmodel on; O model off; and Simodel on;
        if (par.Cmodel == on & par.Omodel == off & par.Simodel == on)
            %
            tmp = sparse(npx+ncx+nsx, npx+ncx+nsx);
            tmp(1:npx, 1:npx) = fxx;
            fxx = tmp;
            % -------------------------------
            % P and C model parameters
            for ju = 1:npx
                for jo = (npx+1):(npx+ncx)
                    kk = kk + 1;
                    fxx(ju, jo) = fxx(ju, jo) + ...
                        cx(idic,ju).'*Wc*cx(idic,jo) + ec.'*Wc*cxx(idic,kk);
                    
                    fxx(jo, ju) = fxx(ju, jo);
                end 
            end
            % -------------------------------
            % C model parameters
            for ju = (npx+1):(npx+ncx)
                for jo = ju:(npx+ncx)
                    kk = kk + 1;
                    fxx(ju, jo) = fxx(ju, jo) + ...
                        cx(idic,ju).'*Wc*cx(idic,jo) + ec.'*Wc*cxx(idic,kk);
                    
                    fxx(jo, ju) = fxx(ju, jo);
                end 
            end
            % -------------------------------
            % P model parameters and Si parameters
            kk = kpp; % starting fro P-P parameters         
            for ju = 1:npx
                for jo = (npx+ncx+1):(npx+ncx+nsx)
                    kk = kk + 1; 
                    fxx(ju, jo) = fxx(ju, jo) + ...
                        sx(isil,ju).'*Ws*sx(isil,jo) + es.'*Ws*sxx(isil,kk);
                    
                    fxx(jo, ju) = fxx(ju, jo);
                end 
            end
            % -----------------------------
            % Only Si parameters
            for ju = (npx+ncx+1):(npx+ncx+nsx)
                for jo = ju:(npx+ncx+nsx)
                    kk = kk + 1;
                    fxx(ju, jo) = fxx(ju, jo) + ...
                        sx(isil,ju).'*Ws*sx(isil,jo) + es.'*Ws*sxx(isil,kk);
                    
                    fxx(jo, ju) = fxx(ju, jo);
                end 
            end  
        end
        % ----------------------------------------------------------------
        % Cmodel on; O model on; and Simodel on;
        if (par.Cmodel == on & par.Omodel == on & par.Simodel == on)
            %
            tmp = sparse(npx+ncx+nox+nsx, npx+ncx+nox+nsx);
            tmp(1:npx, 1:npx) = fxx;
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
                        ox(io2, ju).'*Wo*ox(io2, jo) + eo.'*Wo*oxx(io2, kk);
                    
                    fxx(jo, ju) = fxx(ju, jo);
                end 
            end
            % --------------------------
            % P model parameters and Si parameters
            kk = kpp; % starting fro P-P parameters 
            for ju = 1:npx
                for jo = (npx+ncx+nox+1):(npx+ncx+nox+nsx)
                    kk = kk + 1; 
                    fxx(ju, jo) = fxx(ju, jo) + ...
                        sx(isil,ju).'*Ws*sx(isil,jo) + es.'*Ws*sxx(isil,kk);
                    
                    fxx(jo, ju) = fxx(ju, jo);
                end 
            end
            % -----------------------------
            % Only Si parameters
            for ju = (npx+ncx+nox+1):(npx+ncx+nox+nsx)
                for jo = ju:(npx+ncx+nox+nsx)
                    kk = kk + 1;
                    fxx(ju, jo) = fxx(ju, jo) + ...
                        sx(isil,ju).'*Ws*sx(isil,jo) + es.'*Ws*sxx(isil,kk);
                    
                    fxx(jo, ju) = fxx(ju, jo);
                end 
            end  
        end
    end
end

