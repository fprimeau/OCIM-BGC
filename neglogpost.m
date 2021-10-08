function [f, fx, fxx, data] = neglogpost(x, par)
    global iter
    on = true; off = false;
	current_time = string(datetime('now')) ; %for runtime diagnostic purposes
	fprintf('current time: %s \n',current_time) ;

	nx   = length(x) ; % number of parameters
    dVt  = par.dVt   ;
    M3d  = par.M3d   ;
    iwet = par.iwet  ;
    nwet = par.nwet  ;
    %
    f    = 0 ; % initialize objective function value

    % print and save current parameter values to
    % a file that is used to reset parameters ;
    if iter == 0
        xhat = PrintPara(x, par) ;
    end
    % reset parameters if optimization routine
    % suggests strange parameter values ;
    %if iter < 3
    %    x = ResetPara(x, par) ;
    %end
    % print current parameters
    if iter > 0
        xhat = PrintPara(x, par) ;
    end
    fprintf('current iteration is %d \n',iter) ;
    iter = iter + 1  ;

	%myobjfun = [];
	%myobjfunx = [];

	% do not execute code  if solver suggests very bad values
if iter>1 & iter<5
	ibad = BadStep(x,par);
	% replace resetPara with BadStep(x, par).  use the same stopping criteria, except instead of replacing parameter value, add the pindx to a vector ibad. then if ibad has length > 0, set f to a large number (f=1000); set fx to a vector of zeros the right length; and set fxx to a matrix of correct size. length(fx) = length(x)
	% do this step before printpara, because print para updates the values saved in fxhat
	if ~isempty(ibad)
		load(par.fxhat);
		f = 10000;
		%fx = zeros(length(x), 1);
		%fxx = sparse(nx, nx);
		data = struct;
		fprintf('solver suggested unrealistic parameter values. exiting neglogpost... \n')
		return
	end
end
% save current set of parameter values to a file
if (par.optim == on)
	x0 = x ;
	save(par.fxhat, 'x0','xhat')
	fprintf('saving x0 and xhat to par.fxhat...\n')
end


    %%%%%%%%%%%%%%%%%%   Solve P    %%%%%%%%%%%%%%%%%%%%%%%%
	tic
    idip = find(par.po4raw(iwet) > 0.02) ;
    Wp   = d0(dVt(iwet(idip))/sum(dVt(iwet(idip)))) ;
    mu   = sum(Wp*par.po4raw(iwet(idip)))/sum(diag(Wp)) ;
    var  = sum(Wp*(par.po4raw(iwet(idip))-mu).^2)/sum(diag(Wp)) ;
    Wip  = Wp/var ;

    idop = find(par.dopraw(iwet) > 0.0) ;
    Wp   = d0(dVt(iwet(idop))/sum(dVt(iwet(idop)))) ;
    mu   = sum(Wp*par.dopraw(iwet(idop)))/sum(diag(Wp)) ;
    var  = sum(Wp*(par.dopraw(iwet(idop))-mu).^2)/sum(diag(Wp)) ;
    Wop  = par.pscale*Wp/var ;
    %
    [par, P, Px, Pxx] = eqPcycle(x, par) ;
    DIP = M3d+nan  ;  DIP(iwet) = P(1+0*nwet:1*nwet) ;
    POP = M3d+nan  ;  POP(iwet) = P(1+1*nwet:2*nwet) ;
    DOP = M3d+nan  ;  DOP(iwet) = P(1+2*nwet:3*nwet) ;

    par.Px   = Px  ;
    par.Pxx  = Pxx ;
    par.DIP  = DIP(iwet) ;
    data.DIP = DIP ; data.POP = POP ; data.DOP = DOP ;
    % DIP error
    eip = DIP(iwet(idip)) - par.po4raw(iwet(idip)) ;
    eop = DOP(iwet(idop)) - par.dopraw(iwet(idop)) ;
    f  = f + 0.5*(eip.'*Wip*eip) + 0.5*(eop.'*Wop*eop); %+ prior.bP
	toc
    %%%%%%%%%%%%%%%%%%   End Solve P    %%%%%%%%%%%%%%%%%%%%

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

	%%%%%%%%%%%%%%     Solve for C2P with Cell model   %%%%%%%%%%%%%%%%%%%%%
	if (par.Cellmodel == on)
		tic
		[par, C2P, C2Px, C2Pxx, C2Ppxx] = eqC2Puptake(x, par, data); % this line replaces the rest of this section
		par.C2Px = C2Px;
		par.C2Pxx = C2Pxx;
		par.C2Ppxx = C2Ppxx;

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
		toc
	end

    %%%%%%%%%%%%%%%%%%     Solve C   %%%%%%%%%%%%%%%%%%%%%%%%
    if (par.Cmodel == on)
		tic
        idic = find(par.dicraw(iwet) > 0) ;
        Wic  = d0(dVt(iwet(idic))/sum(dVt(iwet(idic)))) ;
        mu   = sum(Wic*par.dicraw(iwet(idic)))/sum(diag(Wic)) ;
        var  = sum(Wic*(par.dicraw(iwet(idic))-mu).^2)/sum(diag(Wic));
        Wic  = Wic/var  ;

        ialk = find(par.alkraw(iwet)>0) ;
        Wlk  = d0(dVt(iwet(ialk))/sum(dVt(iwet(ialk)))) ;
        mu   = sum(Wlk*par.alkraw(iwet(ialk)))/sum(diag(Wlk)) ;
        var  = sum(Wlk*(par.alkraw(iwet(ialk))-mu).^2)/sum(diag(Wlk));
        Wlk  = Wlk/var  ;

        idoc = find(par.docraw(iwet)>0) ;
        Woc  = d0(dVt(iwet(idoc))/sum(dVt(iwet(idoc)))) ;
        mu   = sum(Woc*par.docraw(iwet(idoc)))/sum(diag(Woc)) ;
        var  = sum(Woc*(par.docraw(iwet(idoc))-mu).^2)/sum(diag(Woc));
        Woc  = par.cscale*Woc/var ;

        [par, C, Cx, Cxx] = eqCcycle(x, par) ;
        DIC = M3d+nan ;  DIC(iwet) = C(0*nwet+1:1*nwet) ;
        POC = M3d+nan ;  POC(iwet) = C(1*nwet+1:2*nwet) ;
        DOC = M3d+nan ;  DOC(iwet) = C(2*nwet+1:3*nwet) ;
        PIC = M3d+nan ;  PIC(iwet) = C(3*nwet+1:4*nwet) ;
        ALK = M3d+nan ;  ALK(iwet) = C(4*nwet+1:5*nwet) ;

        par.DIC = DIC(iwet) ;
        par.DOC = DOC(iwet) ;
        DIC = DIC + par.dicant  ;
        par.Cx   = Cx  ;  par.Cxx  = Cxx ;
        data.DIC = DIC ;  data.POC = POC ;
        data.DOC = DOC ;  data.PIC = PIC ;
        data.ALK = ALK ;
%
%		CNPP = M3d+nan ;  CNPP(iwet) = par.G*par.C2P ;
%		data.CNPP = CNPP;

        % DIC error
        eic = DIC(iwet(idic)) - par.dicraw(iwet(idic)) ;
        eoc = DOC(iwet(idoc)) - par.docraw(iwet(idoc)) ;
        elk = ALK(iwet(ialk)) - par.alkraw(iwet(ialk)) ;
        f   = f + 0.5*(eic.'*Wic*eic) + 0.5*(eoc.'*Woc*eoc) + ...
              0.5*(elk.'*Wlk*elk);
		toc
    end

% debugging gradient --> fixed!
%	[ibadx,ibady] = find(Cx ~= real(Cx));
%	ibadparam = unique(ibady)
%	keyboard;
    %%%%%%%%%%%%%%%%%%   End Solve C    %%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%   Solve O    %%%%%%%%%%%%%%%%%%%%%%%%
    if (par.Omodel == on)
		tic
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
		toc
    end
    %%%%%%%%%%%%%%%%%%   End Solve O    %%%%%%%%%%%%%%%%%%%%
    fprintf('current objective function value is %3.3e \n\n',f)
    if mod(iter, 10) == 0
        save(par.fname, 'data')
    end
    %% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    % calculate gradient
    if (nargout > 1)
        fx = zeros(length(x), 1)   ;
        ipx  = Px(0*nwet+1:nwet,:) ;
        opx  = Px(2*nwet+1:end ,:) ;
        npx = par.npx              ;
        % ---------Pmodel------------------------
        for ji = 1 : npx
            fx(ji) = eip.'*Wip*ipx(idip,ji) + eop.'*Wop*opx(idop,ji); % + 2*10^6(par.bP-1) dobj_fun_dbP = dP_dbP*prior + f*dprior_dbP
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
            ocx = Cx(2*nwet+1 : 3*nwet, :) ;
            lkx = Cx(4*nwet+1 : 5*nwet, :) ;
            ncx = par.ncx   ;
			nbx = par.nbx   ; % treating cell model as C-cycle parameters (replace ncx with ncx+nbx)
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
            for ji = 1 : npx+ncx+nbx+nox	%treating cell model as C-cycle parameters (replace ncx with ncx+nbx)
                fx(ji) = fx(ji) + eo.'*Wo*ox(io2, ji) ;
            end
        end
        % ----------------------------------
        % ----------------------------------
        if (par.Cmodel == on & par.Simodel == on & par.Omodel == off)
            ncx = par.ncx ;
			nbx = par.nbx ; % treating cell model as C-cycle parameters (replace ncx with ncx+nbx)
            nsx = par.nsx ;
            icx = Cx(0*nwet+1 : 1*nwet, :) ;
            ocx = Cx(2*nwet+1 : 3*nwet, :) ;
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
			nbx = par.nbx ; % treating cell model as C-cycle parameters (replace ncx with ncx+nbx)
            nox = par.nox ;
            nsx = par.nsx ;
            icx = Cx(0*nwet+1 : 1*nwet, :) ;
            ocx = Cx(2*nwet+1 : 3*nwet, :) ;
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
            for ji = npx+ncx+nbx+1 : npx+ncx+nbx+nox       % BUG: needs npx+ncx starting indx
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
        ipxx = Pxx(0*nwet+1 : 1*nwet, :) ;
        opxx = Pxx(2*nwet+1 : 3*nwet, :) ;
        % ----------------------------------------------------------------
        kk = 0;
        for ju = 1:npx
            for jo = ju:npx
                kk = kk + 1 ;
                fxx(ju,jo) = ...
                    ipx(idip,ju).'*Wip*ipx(idip,jo) + eip.'*Wip*ipxx(idip,kk) + ...
                    opx(idop,ju).'*Wop*opx(idop,jo) + eop.'*Wop*opxx(idop,kk);
                % Simodel
                if par.Simodel == on
                    sxx = Sixx(1:nwet,:);
                    fxx(ju,jo) = fxx(ju, jo) + ...
                        sx(isil,ju).'*Ws*sx(isil,jo) + es.'*Ws*sxx(isil,kk);
                end
                % Cmodel
                if par.Cmodel == on
                    icxx = Cxx(0*nwet+1 : 1*nwet, :) ;
                    ocxx = Cxx(2*nwet+1 : 3*nwet, :) ;
                    lkxx = Cxx(4*nwet+1 : 5*nwet, :) ;
                    fxx(ju,jo) = fxx(ju, jo) + ...
                        icx(idic,ju).'*Wic*icx(idic,jo) + eic.'*Wic*icxx(idic,kk) + ...
                        ocx(idoc,ju).'*Woc*ocx(idoc,jo) + eoc.'*Woc*ocxx(idoc,kk) + ...
                        lkx(ialk,ju).'*Wlk*lkx(ialk,jo) + elk.'*Wlk*lkxx(ialk,kk) ;
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
			% treating cell model as C-cycle parameters (replace ncx with ncx+nbx)
        if (par.Cmodel == on & par.Omodel == off & par.Simodel == off)
            %
            tmp = sparse(npx+ncx+nbx, npx+ncx+nbx); %
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
            for ju = (npx+1):(npx+ncx)  %to npx+ncx+nbx to include cell params (check order in Cxx)
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

                    fxx(jo,ju) = fxx(ju, jo);
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
			% C and cell model parameters
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
            % C and O model parameters     % need to add cell model with O
            for ju = (npx+1):(npx+ncx+nbx)
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
                for jo = (npx+1):(npx+ncx+nbx)
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
            for ju = (npx+1):(npx+ncx+nbx)
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
                for jo = (npx+1):(npx+ncx+nbx)
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
            for ju = (npx+1):(npx+ncx+nbx)
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
            for ju = (npx+1):(npx+ncx+nbx)
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

		if iter < 5
			load(par.fxhat);
			save(par.fxhat, 'x0','xhat','f','fx','fxx')
		end
	end %(nargout>2)


	% if testing single parameter
	% myobjfun = [myobjfun, f];
	% myobjfunx = [myobjfunx, fx];
	% testname=['objfun_test_' par.fname];
	% save(testname, 'x','objfun','objfunx');

end % end function
