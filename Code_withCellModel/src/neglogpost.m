function [f, fx, fxx, data, xhat] = neglogpost(x, par)
    global iter
    on = true; off = false;
    % print and save current parameter values to
    % a file that is used to reset parameters ;
    if iter == 0
        PrintPar(x, par) ;
    end

    % reset parameters if optimization routine
    % suggests strange parameter values ;
    if iter>0 & iter < 10  %skipping Reset par
        [x, ibad] = ResetPar(x, par) ;
        % Do not execute code if solver suggests very bad values
        % (instead of replacing bad parameter value and solving)
		if ~isempty(ibad)
			load(par.fxhat); % make sure to do this step before saving the reset values in fxhat
			f = 10000;
			fx = xhat.fx;
			fxx = xhat.fxx;
			data = struct;
			fprintf('solver suggested unrealistic parameter values for pindx = ')
            fprintf('%d ... ',ibad)
            fprintf('\n')
            fprintf('Exiting neglogpost... \n')
			return
		end
    end
    % print and save current parameters to a file par.fxhat
    if iter > 0
        PrintPar(x, par) ;    
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
    idip = find(par.po4raw(iwet) > 0.05) ;
    Wp   = d0(dVt(iwet(idip))/sum(dVt(iwet(idip)))) ;
    mu   = sum(Wp*par.po4raw(iwet(idip)))/sum(diag(Wp)) ;
    var  = sum(Wp*(par.po4raw(iwet(idip))-mu).^2)/sum(diag(Wp)) ;
    Wip  = par.dipscale*Wp/var ;

    idop = find(par.dopraw(iwet) > 0.0) ;
    Wp   = d0(dVt(iwet(idop))/sum(dVt(iwet(idop)))) ;
    mu   = sum(Wp*par.dopraw(iwet(idop)))/sum(diag(Wp)) ;
    var  = sum(Wp*(par.dopraw(iwet(idop))-mu).^2)/sum(diag(Wp)) ;
    Wop  = par.dopscale*Wp/var ;
    %
    %tic 
    [par, P, Px, Pxx] = eqPcycle(x, par) ;
    DIP  = M3d+nan  ;  DIP(iwet)  = P(1+0*nwet:1*nwet) ;
    POP  = M3d+nan  ;  POP(iwet)  = P(1+1*nwet:2*nwet) ;
    DOP  = M3d+nan  ;  DOP(iwet)  = P(1+2*nwet:3*nwet) ;
    DOPl = M3d+nan  ;  DOPl(iwet) = P(1+2*nwet:3*nwet) ;
    %toc 
    par.Px   = Px  ;
    par.Pxx  = Pxx ;
    par.DIP  = DIP(iwet) ;
    data.DIP = DIP ; data.POP  = POP  ;
    data.DOP = DOP ; data.DOPl = DOPl ;
    % DIP error
    eip = DIP(iwet(idip)) - par.po4raw(iwet(idip)) ;
    eop = DOP(iwet(idop)) - par.dopraw(iwet(idop)) ;
    f  = f + 0.5*(eip.'*Wip*eip) + 0.5*(eop.'*Wop*eop); 
    
    %%%%%%%%%%%%%%%%%%   End Solve P    %%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%     Solve for C2P with Cell model   %%%%%%%%%%%%%%%%%%%%%
	if (par.Cellmodel == on)
		tic
		[par, C2P, C2Px, C2Pxx, C2Ppxx] = eqC2Puptake(x, par, data);
        par.C2P = C2P ; % vector length nwet
		par.C2Px = C2Px;
		par.C2Pxx = C2Pxx;
		par.C2Ppxx = C2Ppxx;

		%save derivative for sensitivity tests
		data.C2Px = C2Px;

        % save additional cell model fields for analysis
		%data.CellOut = par.CellOut;
		data.CellOut.C2P =  par.CellOut.C2P;
		data.CellOut.N2P = par.CellOut.N2P;
		data.CellOut.C2N = par.CellOut.C2N;
		data.CellOut.LimType = par.CellOut.LimType;
		data.CellOut.r = par.CellOut.r;
		data.CellOut.mu = par.CellOut.mu;
		data.CellOut.E = par.CellOut.E;
		data.CellOut.L = par.CellOut.L;
		data.CellOut.A = par.CellOut.A; % M=A
		data.CellOut.PLip = par.CellOut.PLip;
		data.CellOut.PStor = par.CellOut.PStor;
		data.CellOut.QP = par.CellOut.QP;
		data.CellOut.QC = par.CellOut.QC;
		data.CellOut.dC2P_dDIP = par.CellOut.dC2P_dDIP;

		toc
		fprintf('All of Cell model solved \n ')
	end
       
    %%%%%%%%%%%%%%%%%%   Solve Si       %%%%%%%%%%%%%%%%%%%%
    if (par.Simodel == on)
        isil = find(par.sio4raw(iwet)>0) ;
        Ws   = d0(dVt(iwet(isil))/sum(dVt(iwet(isil)))) ;
        mu   = sum(Ws*par.sDSi(iwet(isil)))/sum(diag(Ws)) ;
        var  = sum(Ws*(par.sDSi(iwet(isil))-mu).^2)/sum(diag(Ws));
        Ws   = Ws/var ;
        %
        [par,Si,Six,Sixx] = eqSicycle(x, par)   ;
        DSi = M3d+nan ;  DSi(iwet) = Si(1:nwet) ;
        bSi = M3d+nan ;  bSi(iwet) = Si(nwet+1:end) ;
        data.DSi = SI ;  data.bSi  = bSi ;  
        % SiO error
        es = DSi(iwet(isil)) - par.sio4raw(iwet(isil)) ;
        f  = f + 0.5*(es.'*Ws*es) ;
    end
    %%%%%%%%%%%%%%%%%%   End Solve Si    %%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%     Solve C   %%%%%%%%%%%%%%%%%%%%%%%%
    if (par.Cmodel == on)
        idic = find(par.dicraw(iwet) > 0) ;
        Wic  = d0(dVt(iwet(idic))/sum(dVt(iwet(idic)))) ;
        mu   = sum(Wic*par.dicraw(iwet(idic)))/sum(diag(Wic)) ;
        var  = sum(Wic*(par.dicraw(iwet(idic))-mu).^2)/sum(diag(Wic));
        Wic  = par.dicscale*Wic/var  ;
        
        ialk = find(par.alkraw(iwet) > 0) ;
        Wlk  = d0(dVt(iwet(ialk))/sum(dVt(iwet(ialk)))) ;
        mu   = sum(Wlk*par.alkraw(iwet(ialk)))/sum(diag(Wlk)) ;
        var  = sum(Wlk*(par.alkraw(iwet(ialk))-mu).^2)/sum(diag(Wlk));
        Wlk  = par.alkscale*Wlk/var  ;
        
        idoc = find(par.docraw(iwet) > 0) ;
        Woc  = d0(dVt(iwet(idoc))/sum(dVt(iwet(idoc)))) ;
        mu   = sum(Woc*par.docraw(iwet(idoc)))/sum(diag(Woc)) ;
        var  = sum(Woc*(par.docraw(iwet(idoc))-mu).^2)/sum(diag(Woc));
        Woc  = par.docscale*Woc/var ;
        %tic 
        [par, C, Cx, Cxx] = eqCcycle_v2(x, par) ;
        DIC  = M3d+nan ;  DIC(iwet)  = C(0*nwet+1:1*nwet) ;
        POC  = M3d+nan ;  POC(iwet)  = C(1*nwet+1:2*nwet) ;
        DOC  = M3d+nan ;  DOC(iwet)  = C(2*nwet+1:3*nwet) ;
        PIC  = M3d+nan ;  PIC(iwet)  = C(3*nwet+1:4*nwet) ;
        ALK  = M3d+nan ;  ALK(iwet)  = C(4*nwet+1:5*nwet) ;
        DOCl = M3d+nan ;  DOCl(iwet) = C(5*nwet+1:6*nwet) ;
        DOCr = M3d+nan ;  DOCr(iwet) = C(6*nwet+1:7*nwet) ;
       % toc

        par.DIC = DIC(iwet) ;
        par.POC = POC(iwet) ;
        par.DOC = DOC(iwet) ;
        par.DOCl = DOCl(iwet) ;
        par.DOCr = DOCr(iwet) ;
        % DIC = DIC + par.dicant  ;
        par.Cx    = Cx   ;  par.Cxx   = Cxx ;
        data.DIC  = DIC  ;  data.POC  = POC ;
        data.DOC  = DOC  ;  data.PIC  = PIC ;
        data.ALK  = ALK  ;  data.DOCr = DOCr ;
        data.DOCl = DOCl ;
        % DIC error
        DOC = DOC + DOCr + DOCl; % sum of labile and refractory DOC ;
        eic = DIC(iwet(idic)) - par.dicraw(iwet(idic)) ;
        eoc = DOC(iwet(idoc)) - par.docraw(iwet(idoc)) ;
        elk = ALK(iwet(ialk)) - par.alkraw(iwet(ialk)) ;
        f   = f + 0.5*(eic.'*Wic*eic) + 0.5*(eoc.'*Woc*eoc) + ...
              0.5*(elk.'*Wlk*elk);
    end
    %%%%%%%%%%%%%%%%%%   End Solve C    %%%%%%%%%%%%%%%%%%%
      
    %%%%%%%%%%%%%%%%%%   Solve O    %%%%%%%%%%%%%%%%%%%%%%%%
    if (par.Omodel == on)
        io2 = find(par.o2raw(iwet)>0) ;
        Wo  = d0(dVt(iwet(io2))/sum(dVt(iwet(io2)))) ;
        mu  = sum(Wo*par.o2raw(iwet(io2)))/sum(diag(Wo)) ;
        var = sum(Wo*(par.o2raw(iwet(io2))-mu).^2)/sum(diag(Wo)) ;
        Wo  = par.o2scale*Wo/var ;
        %
       % tic 
        [par, O, Ox, Oxx] = eqOcycle_v2(x, par) ;
       % toc 
        O2 = M3d+nan ;  O2(iwet) = O ;
        data.O2 = O2 ;
        eo = O2(iwet(io2)) - par.o2raw(iwet(io2)) ;
        f  = f + 0.5*(eo.'*Wo*eo)   ;
    end
    %%%%%%%%%%%%%%%%%%   End Solve O    %%%%%%%%%%%%%%%%%%%%
    % --------- Save all parameter values inlcuding fixed values in xhat ------
	% this is to make sure all parameter values are consistent when rerunning model
    % parameter_names_P = {'sigP', 'Q10P', 'kdP','bP_T', 'bP','alpha','beta'};
    % parameter_names_Cell = {'Q10Photo','fRibE','fStorage','kST0','PStor_rCutoff','PStor_scale','PLip_PCutoff','PLip_scale','alphaS','gammaDNA'};
    % parameter_names_C = {'sigC', 'kru', 'krd', 'etau', 'etad', 'bC_T', 'bC', 'd', 'Q10C', 'kdC', 'R_Si', 'rR', 'cc', 'dd'};
    % parameter_names_O = {'O2C_T', 'rO2C'};
    % optimizable_parameter_list = [parameter_names_P, parameter_names_Cell, parameter_names_C, parameter_names_O];

    % fixed parameters
	allparams.kPIC      = par.kPIC      ;   % PIC dissolution timescale [1/s]
    allparams.kappa_p 	= par.kappa_p   ;
    allparams.kappa_l 	= par.kappa_l   ;   % labile DOM remin timescale [1/s]
	allparams.kappa_g 	= par.kappa_g   ;
    allparams.gamma     = par.gamma     ;
    % P model parameters
	allparams.sigP 	    = par.sigP     ;
	allparams.Q10P		= par.Q10P      ;
	allparams.kdP 		= par.kdP       ;
	allparams.bP_T 		= par.bP_T      ;
	allparams.bP 		= par.bP        ;
	allparams.alpha 	= par.alpha    ;
	allparams.beta 		= par.beta      ;
	if (par.Cmodel == on)
        parameter_names_C = {'sigC', 'kru', 'krd', 'etau', 'etad', 'bC_T', 'bC', 'd', 'Q10C', 'kdC', 'R_Si', 'rR'}; %excluding cc and dd
        for ii = 1:length(parameter_names_C)
            allparams.(parameter_names_C{ii}) = par.(parameter_names_C{ii});
        end
        %C2P model params
        if par.C2P_Tzmodel
            allparams.ccT = par.ccT; 
            allparams.ddT = par.ddT;
        else
            allparams.cc = par.cc; 
            allparams.dd = par.dd;
        end
	end
	if (par.Omodel == on)
		allparams.O2C_T = par.O2C_T ;
		allparams.rO2C = par.rO2C  ;
	end
	if (par.Simodel==on)
		allparams.dsi =  par.dsi   ;
		allparams.at = par.at    ;
		allparams.bt = par.bt    ;
		allparams.iaa = par.iaa   ;
		allparams.bb = par.bb  ;
	end
	if (par.Cellmodel==on)
		%pull params from BIO substructure
		%allparams.BIO = par.BIO;
		BIOfnames = fieldnames(par.BIO);
		for ii = 1:length(BIOfnames)
			allparams.(BIOfnames{ii}) = par.BIO.(BIOfnames{ii});
		end
		clear BIOfnames
	end
	xhat.allparams = allparams;
    %--------------------------------------------------------
    fprintf('current objective function value is %3.3e \n\n',f) 

    % Save output every (1) iterations during optimization
    if mod(iter, 1) == 0
        fprintf('saving data to %s  ...\n', par.fname) 
        save(par.fname, 'data')
        x0 = real(x) ;   
        fprintf('saving xhat to %s  ...\n', par.fxhat)
        save(par.fxhat, 'xhat','x0')  % save x0 for ResetPar; this xhat also includes allparams (allparams not saved in PrintPar)
        clear x0;
		
	    % save(par.fxpar, 'par' , '-v7.3')
    end 
    %% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    % calculate gradient
    if (nargout > 1)
        fx = zeros(length(x), 1)   ;
        npx = par.npx  ;
        if npx > 0
            ipx  = Px(0*nwet+1:nwet,:) ;
            opx  = Px(2*nwet+1:end ,:) ;
            % ---------------------------------
            for ji = 1 : npx
                fx(ji) = eip.'*Wip*ipx(idip,ji) + eop.'*Wop*opx(idop,ji);
            end
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
            icx = Cx(0*nwet+1 : 1*nwet, :) ;
            ocx = Cx(2*nwet+1 : 3*nwet, :) + Cx(5*nwet+1 : 6*nwet, :) + Cx(6*nwet+1 : 7*nwet, :);
            lkx = Cx(4*nwet+1 : 5*nwet, :) ;
            ncx = par.ncx       ;
            nbx = par.nbx       ;
            % --------------------------
            for ji = 1 : npx+ncx+nbx
                fx(ji) = eic.'*Wic*icx(idic, ji) + ...
                         eoc.'*Woc*ocx(idoc, ji) + ...
                         elk.'*Wlk*lkx(ialk, ji) + ...
                         fx(ji);
            end
        end 
        % ----------------------------------
        % ----------------------------------
        if (par.Omodel == on & par.Simodel == off)
            ox  = Ox(1:nwet, :) ;
            nox = par.nox       ;
            % --------------------------
            for ji = 1 : npx+ncx+nbx+nox
                fx(ji) = fx(ji) + eo.'*Wo*ox(io2, ji) ;
            end
        end 
        % ----------------------------------
        % ----------------------------------
        if (par.Cmodel == on & par.Simodel == on & par.Omodel == off)
            ncx = par.ncx ;
            nbx = par.nbx ;
            nsx = par.nsx ;
            icx = Cx(0*nwet+1 : 1*nwet, :) ;
            ocx = Cx(2*nwet+1 : 3*nwet, :) + Cx(5*nwet+1 : 6*nwet, :) + Cx(6*nwet+1 : 7*nwet, :);
            lkx = Cx(4*nwet+1 : 5*nwet, :) ;
            sx  = Six(1:nwet, :) ;
            % --------------------------
            for ji = 1 : npx
                fx(ji) = eic.'*Wic*icx(idic, ji) + ...
                         eoc.'*Woc*ocx(idoc, ji) + ...
                         elk.'*Wlk*lkx(ialk, ji) + ...
                         es.'*Ws*sx(isil, ji) + fx(ji) ;
            end
            % --------------------------
            for ji = npx+1 : npx+ncx+nbx
                fx(ji) = eic.'*Wic*icx(idic, ji) + ...
                         eoc.'*Woc*ocx(idoc, ji) + ...
                         elk.'*Wlk*lkx(ialk, ji) + ...
                         fx(ji); 
            end
            % --------------------------
            for ji = npx+ncx+nbx+1 : npx+ncx+nbx+nsx
                fx(ji) = fx(ji) + es.'*Ws*sx(isil, ji) ;
            end
        end 
        % ----------------------------------
        % ----------------------------------
        if (par.Cmodel == on & par.Simodel == on & par.Omodel == on)
            ncx = par.ncx ;
            nbx = par.nbx ;
            nox = par.nox ;
            nsx = par.nsx ;
            icx = Cx(0*nwet+1 : 1*nwet, :) ;
            ocx = Cx(2*nwet+1 : 3*nwet, :) + Cx(5*nwet+1 : 6*nwet, :) + Cx(6*nwet+1 : 7*nwet, :);
            lkx = Cx(4*nwet+1 : 5*nwet, :) ;
            ox  = Ox(1:nwet,:) ;
            sx  = Six(1:nwet,:);
            % --------------------------
            for ji = 1 : npx
                fx(ji) = eic.'*Wic*icx(idic, ji) + ...
                         eoc.'*Woc*ocx(idoc, ji) + ...
                         elk.'*Wlk*lkx(ialk, ji) + ...
                         eo.'*Wo*ox(io2 , ji) + ...
                         es.'*Ws*sx(isil, ji) + fx(ji) ;
            end
            % --------------------------
            for ji = npx+1 : npx+ncx+nbx
                fx(ji) = eic.'*Wic*icx(idic, ji) + ...
                         eoc.'*Woc*ocx(idoc, ji) + ...
                         elk.'*Wlk*lkx(ialk, ji) + ...
                         fx(ji);
            end
            % --------------------------
            for ji = npx+1 : npx+ncx+nbx+nox
                fx(ji) = fx(ji) + eo.'*Wo*ox(io2, ji) ;
            end
            % --------------------------
            for ji = npx+ncx+nbx+nox+1 : npx+ncx+nbx+nox+nsx
                fx(ji) = fx(ji) + es.'*Ws*sx(isil, ji) ;
            end
        end 
        % ----------------------------------
    end 
    %
    if (nargout>2)
        fxx = sparse(npx, npx)  ;
        if par.npx > 0
            ipxx = Pxx(0*nwet+1 : 1*nwet, :) ;
            opxx = Pxx(2*nwet+1 : 3*nwet, :) ;
        end
        if par.Cmodel == on & (par.npx + par.ncx + par.nbx > 0)
			icxx = Cxx(0*nwet+1 : 1*nwet, :) ;
            ocxx = Cxx(2*nwet+1 : 3*nwet, :) + Cxx(5*nwet+1 : 6*nwet, :) + Cxx(6*nwet+1 : 7*nwet, :);
            lkxx = Cxx(4*nwet+1 : 5*nwet, :) ; 
		end
		if par.Omodel == on & (par.npx + par.ncx + par.nbx + par.nox > 0)
			oxx = Oxx(1:nwet,:);
		end
		if par.Simodel == on & (par.npx + par.nsx > 0)
			sxx = Sixx(1:nwet,:);
		end
        % ----------------------------------------------------------------
        % P model parameters
        kk = 0;
        for ju = 1:npx
            for jo = ju:npx
                kk = kk + 1 ;
                fxx(ju,jo) = ...
                    ipx(idip,ju).'*Wip*ipx(idip,jo) + eip.'*Wip*ipxx(idip,kk) + ...
                    opx(idop,ju).'*Wop*opx(idop,jo) + eop.'*Wop*opxx(idop,kk);
                % Simodel
                if (par.Simodel == on) & (par.npx + par.nsx > 0)
                    % sxx = Sixx(1:nwet,:);            
                    fxx(ju,jo) = fxx(ju, jo) + ...
                        sx(isil,ju).'*Ws*sx(isil,jo) + es.'*Ws*sxx(isil,kk);
                end 
                % Cmodel
                if (par.Cmodel == on) & (par.npx + par.ncx + par.nbx > 0)
                    % icxx = Cxx(0*nwet+1 : 1*nwet, :) ;
                    % ocxx = Cxx(2*nwet+1 : 3*nwet, :) + Cxx(5*nwet+1 : 6*nwet, :) + Cxx(6*nwet+1 : 7*nwet, :);
                    % lkxx = Cxx(4*nwet+1 : 5*nwet, :) ; 
                    fxx(ju,jo) = fxx(ju, jo) + ...
                        icx(idic,ju).'*Wic*icx(idic,jo) + eic.'*Wic*icxx(idic,kk) + ...
                        ocx(idoc,ju).'*Woc*ocx(idoc,jo) + eoc.'*Woc*ocxx(idoc,kk) + ...
                        lkx(ialk,ju).'*Wlk*lkx(ialk,jo) + elk.'*Wlk*lkxx(ialk,kk) ;
                end
                % Omodel
                if (par.Omodel == on) & (par.npx + par.ncx + par.nbx + par.nox > 0)
                    % oxx = Oxx(1:nwet,:);            
                    fxx(ju,jo) = fxx(ju, jo) + ...
                        ox(io2,ju).'*Wo*ox(io2,jo) + eo.'*Wo*oxx(io2,kk);
                end 
                % make Hessian symetric;
                fxx(jo, ju) = fxx(ju, jo);  
            end 
        end
        kpp = kk;
        % ----------------------------------------------------------------
        % C model parameters; NOTE: treats cell model as C-cycle parameters (replace ncx with ncx+nbx)
        if (par.Cmodel == on & par.Omodel == off & par.Simodel == off)
            %
            tmp = sparse(npx+ncx+nbx, npx+ncx+nbx);
            tmp(1:npx,1:npx) = fxx;
            fxx = tmp;
            % -------------------------------
            % P model parameters and C parameters
            for ju = 1:npx
                for jo = (npx+1):(npx+ncx)
                    kk = kk + 1;
                    fxx(ju, jo) = fxx(ju, jo) + ...
                        icx(idic,ju).'*Wic*icx(idic,jo) + eic.'*Wic*icxx(idic,kk) + ...
                        ocx(idoc,ju).'*Woc*ocx(idoc,jo) + eoc.'*Woc*ocxx(idoc,kk) + ...
                        lkx(ialk,ju).'*Wlk*lkx(ialk,jo) + elk.'*Wlk*lkxx(ialk,kk) ;

                    fxx(jo,ju) = fxx(ju, jo);
                end 
            end
            % P model parameters and Cell parameters
            for ju = 1:npx
                for jo = (npx+ncx+1):(npx+ncx+nbx)
                    kk = kk + 1;
                    fxx(ju, jo) = fxx(ju, jo) + ...
                        icx(idic,ju).'*Wic*icx(idic,jo) + eic.'*Wic*icxx(idic,kk) + ...
                        ocx(idoc,ju).'*Woc*ocx(idoc,jo) + eoc.'*Woc*ocxx(idoc,kk) + ...
                        lkx(ialk,ju).'*Wlk*lkx(ialk,jo) + elk.'*Wlk*lkxx(ialk,kk) ;

                    fxx(jo,ju) = fxx(ju, jo);
                end 
            end
            % -------------------------------
            % Only C parameters
            for ju = (npx+1):(npx+ncx)
                for jo = ju:(npx+ncx)
                    kk = kk + 1;
                    fxx(ju, jo) = fxx(ju, jo) + ...
                        icx(idic,ju).'*Wic*icx(idic,jo) + eic.'*Wic*icxx(idic,kk) + ...
                        ocx(idoc,ju).'*Woc*ocx(idoc,jo) + eoc.'*Woc*ocxx(idoc,kk) + ...
                        lkx(ialk,ju).'*Wlk*lkx(ialk,jo) + elk.'*Wlk*lkxx(ialk,kk) ;

                    fxx(jo, ju) = fxx(ju, jo);
                end 
            end
            % -------------------------------
            % C model parameters and Cell model parameters
            for ju = (npx+1):(npx+ncx)
                for jo = (npx+ncx+1):(npx+ncx+nbx)
                    kk = kk + 1;
                    fxx(ju, jo) = fxx(ju, jo) + ...
                        icx(idic,ju).'*Wic*icx(idic,jo) + eic.'*Wic*icxx(idic,kk) + ...
                        ocx(idoc,ju).'*Woc*ocx(idoc,jo) + eoc.'*Woc*ocxx(idoc,kk) + ...
                        lkx(ialk,ju).'*Wlk*lkx(ialk,jo) + elk.'*Wlk*lkxx(ialk,kk) ;

                    fxx(jo, ju) = fxx(ju, jo);
                end 
            end
            % -------------------------------
            % Only Cell parameters
            for ju = (npx+ncx+1):(npx+ncx+nbx)
                for jo = ju:(npx+ncx+nbx)
                    kk = kk + 1;
                    fxx(ju, jo) = fxx(ju, jo) + ...
                        icx(idic,ju).'*Wic*icx(idic,jo) + eic.'*Wic*icxx(idic,kk) + ...
                        ocx(idoc,ju).'*Woc*ocx(idoc,jo) + eoc.'*Woc*ocxx(idoc,kk) + ...
                        lkx(ialk,ju).'*Wlk*lkx(ialk,jo) + elk.'*Wlk*lkxx(ialk,kk) ;

                    fxx(jo, ju) = fxx(ju, jo);
                end 
            end
        end
        % ----------------------------------------------------------------
        % C and O model
        if (par.Cmodel == on & par.Omodel == on & par.Simodel == off)
            %
            tmp = sparse(npx+ncx+nbx+nox, npx+ncx+nbx+nox);
            tmp(1:npx,1:npx) = fxx;
            fxx = tmp;
            % -------------------------------
            % P and C model parameters
            for ju = 1:npx
                for jo = (npx+1):(npx+ncx)
                    kk = kk + 1;
                    fxx(ju, jo) = ...
                        fxx(ju, jo) + ...
                        icx(idic,ju).'*Wic*icx(idic,jo) + eic.'*Wic*icxx(idic,kk) + ...
                        ocx(idoc,ju).'*Woc*ocx(idoc,jo) + eoc.'*Woc*ocxx(idoc,kk) + ...
                        lkx(ialk,ju).'*Wlk*lkx(ialk,jo) + elk.'*Wlk*lkxx(ialk,kk) + ...
                        ox(io2,ju).'*Wo*ox(io2,jo) + eo.'*Wo*oxx(io2,kk);

                    fxx(jo, ju) = fxx(ju, jo);
                end 
            end
            % -------------------------------
            % P and Cell model parameters
            for ju = 1:npx
                for jo = (npx+ncx+1):(npx+ncx+nbx)
                    kk = kk + 1;
                    fxx(ju, jo) = ...
                        fxx(ju, jo) + ...
                        icx(idic,ju).'*Wic*icx(idic,jo) + eic.'*Wic*icxx(idic,kk) + ...
                        ocx(idoc,ju).'*Woc*ocx(idoc,jo) + eoc.'*Woc*ocxx(idoc,kk) + ...
                        lkx(ialk,ju).'*Wlk*lkx(ialk,jo) + elk.'*Wlk*lkxx(ialk,kk) + ...
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
                        icx(idic,ju).'*Wic*icx(idic,jo) + eic.'*Wic*icxx(idic,kk) + ...
                        ocx(idoc,ju).'*Woc*ocx(idoc,jo) + eoc.'*Woc*ocxx(idoc,kk) + ...
                        lkx(ialk,ju).'*Wlk*lkx(ialk,jo) + elk.'*Wlk*lkxx(ialk,kk) + ...
                        ox(io2,ju).'*Wo*ox(io2,jo) + eo.'*Wo*oxx(io2,kk);
                    
                    fxx(jo, ju) = fxx(ju, jo);
                end 
            end
            % C and Cell model parameters
            for ju = (npx+1):(npx+ncx)
                for jo = (npx+ncx+1):(npx+ncx+nbx)
                    kk = kk + 1;
                    fxx(ju, jo) = ...
                        fxx(ju, jo) + ...
                        icx(idic,ju).'*Wic*icx(idic,jo) + eic.'*Wic*icxx(idic,kk) + ...
                        ocx(idoc,ju).'*Woc*ocx(idoc,jo) + eoc.'*Woc*ocxx(idoc,kk) + ...
                        lkx(ialk,ju).'*Wlk*lkx(ialk,jo) + elk.'*Wlk*lkxx(ialk,kk) + ...
                        ox(io2,ju).'*Wo*ox(io2,jo) + eo.'*Wo*oxx(io2,kk);
                    
                    fxx(jo, ju) = fxx(ju, jo);
                end 
            end
            % -------------------------------
            % Only Cell parameters
            for ju = (npx+ncx+1):(npx+ncx+nbx)
                for jo = ju:(npx+ncx+nbx)
                    kk = kk + 1;
                    fxx(ju, jo) = ...
                        fxx(ju, jo) + ...
                        icx(idic,ju).'*Wic*icx(idic,jo) + eic.'*Wic*icxx(idic,kk) + ...
                        ocx(idoc,ju).'*Woc*ocx(idoc,jo) + eoc.'*Woc*ocxx(idoc,kk) + ...
                        lkx(ialk,ju).'*Wlk*lkx(ialk,jo) + elk.'*Wlk*lkxx(ialk,kk) + ...
                        ox(io2,ju).'*Wo*ox(io2,jo) + eo.'*Wo*oxx(io2,kk);
                    
                    fxx(jo, ju) = fxx(ju, jo);
                end 
            end
            % -------------------------------
            % P and O model parameters
            for ju = 1:npx
                for jo = (npx+ncx+nbx+1):(npx+ncx+nbx+nox)
                    kk = kk + 1;
                    fxx(ju, jo) = fxx(ju, jo) + ...
                        ox(io2,ju).'*Wo*ox(io2,jo) + eo.'*Wo*oxx(io2,kk);
                    
                    fxx(jo, ju) = fxx(ju, jo);
                end 
            end
            % -------------------------------
            % C and O model parameters
            for ju = (npx+1):(npx+ncx)
                for jo = (npx+ncx+nbx+1):(npx+ncx+nbx+nox)
                    kk = kk + 1;
                    fxx(ju, jo) = fxx(ju, jo) + ...
                        ox(io2,ju).'*Wo*ox(io2,jo) + eo.'*Wo*oxx(io2,kk);
                    
                    fxx(jo, ju) = fxx(ju, jo);
                end 
            end
            % -------------------------------
            % Cell and O model parameters
            for ju = (npx+ncx+1):(npx+ncx+nbx)
                for jo = (npx+ncx+nbx+1):(npx+ncx+nbx+nox)
                    kk = kk + 1;
                    fxx(ju, jo) = fxx(ju, jo) + ...
                        ox(io2,ju).'*Wo*ox(io2,jo) + eo.'*Wo*oxx(io2,kk);
                    
                    fxx(jo, ju) = fxx(ju, jo);
                end 
            end
            % -------------------------------
            % O model parameters
            for ju = (npx+ncx+nbx+1):(npx+ncx+nbx+nox)
                for jo = ju:(npx+ncx+nbx+nox)
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
            tmp = sparse(npx+ncx+nbx+nsx, npx+ncx+nbx+nsx);
            tmp(1:npx, 1:npx) = fxx;
            fxx = tmp;
            % -------------------------------
            % P and C model parameters
            for ju = 1:npx
                for jo = (npx+1):(npx+ncx)
                    kk = kk + 1;
                    fxx(ju, jo) = fxx(ju, jo) + ...
                        icx(idic,ju).'*Wic*icx(idic,jo) + eic.'*Wic*icxx(idic,kk) + ...
                        ocx(idoc,ju).'*Woc*ocx(idoc,jo) + eoc.'*Woc*ocxx(idoc,kk) + ...
                        lkx(ialk,ju).'*Wlk*lkx(ialk,jo) + elk.'*Wlk*lkxx(ialk,kk) ;

                    fxx(jo, ju) = fxx(ju, jo);
                end 
            end
            % P and Cell model parameters
            for ju = 1:npx
                for jo = (npx+ncx+1):(npx+ncx+nbx)
                    kk = kk + 1;
                    fxx(ju, jo) = fxx(ju, jo) + ...
                        icx(idic,ju).'*Wic*icx(idic,jo) + eic.'*Wic*icxx(idic,kk) + ...
                        ocx(idoc,ju).'*Woc*ocx(idoc,jo) + eoc.'*Woc*ocxx(idoc,kk) + ...
                        lkx(ialk,ju).'*Wlk*lkx(ialk,jo) + elk.'*Wlk*lkxx(ialk,kk) ;

                    fxx(jo, ju) = fxx(ju, jo);
                end 
            end
            % -------------------------------
            % C model parameters
            for ju = (npx+1):(npx+ncx)
                for jo = ju:(npx+ncx)
                    kk = kk + 1;
                    fxx(ju, jo) = fxx(ju, jo) + ...
                        icx(idic,ju).'*Wic*icx(idic,jo) + eic.'*Wic*icxx(idic,kk) + ...
                        ocx(idoc,ju).'*Woc*ocx(idoc,jo) + eoc.'*Woc*ocxx(idoc,kk) + ...
                        lkx(ialk,ju).'*Wlk*lkx(ialk,jo) + elk.'*Wlk*lkxx(ialk,kk);
                    
                    fxx(jo, ju) = fxx(ju, jo);
                end 
            end
            % C and Cell model parameters
            for ju = (npx+1):(npx+ncx)
                for jo = (npx+ncx+1):(npx+ncx+nbx)
                    kk = kk + 1;
                    fxx(ju, jo) = fxx(ju, jo) + ...
                        icx(idic,ju).'*Wic*icx(idic,jo) + eic.'*Wic*icxx(idic,kk) + ...
                        ocx(idoc,ju).'*Woc*ocx(idoc,jo) + eoc.'*Woc*ocxx(idoc,kk) + ...
                        lkx(ialk,ju).'*Wlk*lkx(ialk,jo) + elk.'*Wlk*lkxx(ialk,kk);
                    
                    fxx(jo, ju) = fxx(ju, jo);
                end 
            end
            % Cell model parameters
            for ju = (npx+ncx+1):(npx+ncx+nbx)
                for jo = ju:(npx+ncx+nbx)
                    kk = kk + 1;
                    fxx(ju, jo) = fxx(ju, jo) + ...
                        icx(idic,ju).'*Wic*icx(idic,jo) + eic.'*Wic*icxx(idic,kk) + ...
                        ocx(idoc,ju).'*Woc*ocx(idoc,jo) + eoc.'*Woc*ocxx(idoc,kk) + ...
                        lkx(ialk,ju).'*Wlk*lkx(ialk,jo) + elk.'*Wlk*lkxx(ialk,kk);
                    
                    fxx(jo, ju) = fxx(ju, jo);
                end 
            end
            % -------------------------------
            % P model parameters and Si parameters
            kk = kpp; % starting fro P-P parameters         
            for ju = 1:npx
                for jo = (npx+ncx+nbx+1):(npx+ncx+nbx+nsx)
                    kk = kk + 1; 
                    fxx(ju, jo) = fxx(ju, jo) + ...
                        sx(isil,ju).'*Ws*sx(isil,jo) + es.'*Ws*sxx(isil,kk);
                    
                    fxx(jo, ju) = fxx(ju, jo);
                end 
            end
            % -----------------------------
            % Only Si parameters
            for ju = (npx+ncx+nbx+1):(npx+ncx+nbx+nsx)
                for jo = ju:(npx+ncx+nbx+nsx)
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
            tmp = sparse(npx+ncx+nbx+nox+nsx, npx+ncx+nbx+nox+nsx);
            tmp(1:npx, 1:npx) = fxx;
            fxx = tmp;
            % -------------------------------
            % P and C model parameters
            for ju = 1:npx
                for jo = (npx+1):(npx+ncx)
                    kk = kk + 1;
                    fxx(ju, jo) = ...
                        fxx(ju, jo) + ...
                        icx(idic,ju).'*Wic*icx(idic,jo) + eic.'*Wic*icxx(idic,kk) + ...
                        ocx(idoc,ju).'*Woc*ocx(idoc,jo) + eoc.'*Woc*ocxx(idoc,kk) + ... 
                        lkx(ialk,ju).'*Wlk*lkx(ialk,jo) + elk.'*Wlk*lkxx(ialk,kk) + ...
                        ox(io2,ju).'*Wo*ox(io2,jo) + eo.'*Wo*oxx(io2,kk);
                    
                    fxx(jo, ju) = fxx(ju, jo);
                end 
            end
            % P and Cell model parameters
            for ju = 1:npx
                for jo = (npx+ncx+1):(npx+ncx+nbx)
                    kk = kk + 1;
                    fxx(ju, jo) = ...
                        fxx(ju, jo) + ...
                        icx(idic,ju).'*Wic*icx(idic,jo) + eic.'*Wic*icxx(idic,kk) + ...
                        ocx(idoc,ju).'*Woc*ocx(idoc,jo) + eoc.'*Woc*ocxx(idoc,kk) + ... 
                        lkx(ialk,ju).'*Wlk*lkx(ialk,jo) + elk.'*Wlk*lkxx(ialk,kk) + ...
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
                        icx(idic,ju).'*Wic*icx(idic,jo) + eic.'*Wic*icxx(idic,kk) + ...
                        ocx(idoc,ju).'*Woc*ocx(idoc,jo) + eoc.'*Woc*ocxx(idoc,kk) + ... 
                        lkx(ialk,ju).'*Wlk*lkx(ialk,jo) + elk.'*Wlk*lkxx(ialk,kk) + ...
                        ox(io2,ju).'*Wo*ox(io2,jo) + eo.'*Wo*oxx(io2,kk);
                    
                    fxx(jo, ju) = fxx(ju, jo);
                end 
            end
            % C and Cell model parameters
            for ju = (npx+1):(npx+ncx)
                for jo = (npx+ncx+1):(npx+ncx+nbx)
                    kk = kk + 1;
                    fxx(ju, jo) = ...
                        fxx(ju, jo) + ...
                        icx(idic,ju).'*Wic*icx(idic,jo) + eic.'*Wic*icxx(idic,kk) + ...
                        ocx(idoc,ju).'*Woc*ocx(idoc,jo) + eoc.'*Woc*ocxx(idoc,kk) + ... 
                        lkx(ialk,ju).'*Wlk*lkx(ialk,jo) + elk.'*Wlk*lkxx(ialk,kk) + ...
                        ox(io2,ju).'*Wo*ox(io2,jo) + eo.'*Wo*oxx(io2,kk);
                    
                    fxx(jo, ju) = fxx(ju, jo);
                end 
            end
            % Cell model parameters
            for ju = (npx+ncx+1):(npx+ncx+nbx)
                for jo = ju:(npx+ncx+nbx)
                    kk = kk + 1;
                    fxx(ju, jo) = ...
                        fxx(ju, jo) + ...
                        icx(idic,ju).'*Wic*icx(idic,jo) + eic.'*Wic*icxx(idic,kk) + ...
                        ocx(idoc,ju).'*Woc*ocx(idoc,jo) + eoc.'*Woc*ocxx(idoc,kk) + ... 
                        lkx(ialk,ju).'*Wlk*lkx(ialk,jo) + elk.'*Wlk*lkxx(ialk,kk) + ...
                        ox(io2,ju).'*Wo*ox(io2,jo) + eo.'*Wo*oxx(io2,kk);
                    
                    fxx(jo, ju) = fxx(ju, jo);
                end 
            end
            % -------------------------------
            % P and O model parameters
            for ju = 1:npx
                for jo = (npx+ncx+nbx+1):(npx+ncx+nbx+nox)
                    kk = kk + 1;
                    fxx(ju, jo) = fxx(ju, jo) + ...
                        ox(io2,ju).'*Wo*ox(io2,jo) + eo.'*Wo*oxx(io2,kk);
                    
                    fxx(jo, ju) = fxx(ju, jo);
                end 
            end
            % -------------------------------
            % C and O model parameters
            for ju = (npx+1):(npx+ncx)
                for jo = (npx+ncx+nbx+1):(npx+ncx+nbx+nox)
                    kk = kk + 1;
                    fxx(ju, jo) = fxx(ju, jo) + ...
                        ox(io2,ju).'*Wo*ox(io2,jo) + eo.'*Wo*oxx(io2,kk);
                    
                    fxx(jo, ju) = fxx(ju, jo);
                end 
            end
            % Cell and O model parameters
            for ju = (npx+ncx+1):(npx+ncx+nbx)
                for jo = (npx+ncx+nbx+1):(npx+ncx+nbx+nox)
                    kk = kk + 1;
                    fxx(ju, jo) = fxx(ju, jo) + ...
                        ox(io2,ju).'*Wo*ox(io2,jo) + eo.'*Wo*oxx(io2,kk);
                    
                    fxx(jo, ju) = fxx(ju, jo);
                end 
            end
            % -------------------------------
            % O model parameters
            for ju = (npx+ncxnbx+1):(npx+ncx+nbx+nox)
                for jo = ju:(npx+ncx+nbx+nox)
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
                for jo = (npx+ncx+nbx+nox+1):(npx+ncx+nbx+nox+nsx)
                    kk = kk + 1; 
                    fxx(ju, jo) = fxx(ju, jo) + ...
                        sx(isil,ju).'*Ws*sx(isil,jo) + es.'*Ws*sxx(isil,kk);
                    
                    fxx(jo, ju) = fxx(ju, jo);
                end 
            end
            % -----------------------------
            % Only Si parameters
            for ju = (npx+ncx+nbx+nox+1):(npx+ncx+nbx+nox+nsx)
                for jo = ju:(npx+ncx+nbx+nox+nsx)
                    kk = kk + 1;
                    fxx(ju, jo) = fxx(ju, jo) + ...
                        sx(isil,ju).'*Ws*sx(isil,jo) + es.'*Ws*sxx(isil,kk);
                    
                    fxx(jo, ju) = fxx(ju, jo);
                end 
            end  
        end
    end
end

