function [par,Si,Six,Sixx] = eqSicycle(x, par)
    on = true; off = false;
    par.nx = length(x);
    pindx = par.pindx;
    M3d   =  par.M3d;
    iwet  = par.iwet;
    nwet  = par.nwet;
    I     = speye(nwet); % make an identity matrix;
                         % unpack the parameters to be optimized
    % dsi
    if (par.opt_dsi == on)
        ldsi = x(pindx.ldsi);
        par.dsi = exp(ldsi);
    end
    % at
    if (par.opt_at == on)
        lat = x(pindx.lat);
        par.at = exp(lat);
    end
    % bt
    if (par.opt_bt == on)
        lbt = x(pindx.lbt);
        par.bt = exp(lbt);
    end
    % aa
    if (par.opt_aa == on)
        par.aa = x(pindx.aa);
    end
    % bb
    if (par.opt_bb == on)
        lbb = x(pindx.lbb);
        par.bb = exp(lbb);
    end
    %
    % ++++++++++++++++++++++++++++++++++++++++++++++++++
    at  = par.at;
    bt  = par.bt;
    aa  = par.aa;
    bb  = par.bb;
    dsi = par.dsi;
    kappa_gs = par.kappa_gs;
    T   = par.modT(iwet) + 273.15;
    kappa_si = at * exp(-bt./T);
    par.kappa_si = kappa_si;
    % ++++++++++++++++++++++++++++++++++++++++++++++++++
    %
    % fixed parameters
    SI4 = par.SIL(iwet);
    SILbar = par.SILbar*M3d(iwet);
    %
    %%%%%%%%%%%%%%%%%%% create Si2P %%%%%%%%%%%%%%%%%%%%%
    smsk = M3d;
    smsk(:,:,3:end) = 0;
    isrf = find(smsk(iwet));
    dVs = par.dVt(iwet(isrf));
    surface_mean = @(x) sum(x(isrf).*dVs)/sum(dVs);
    
    % % commentmpute the mean of the regressor variable
    Z = SI4;
    mu = surface_mean(Z);
    Delta = sqrt(surface_mean((Z-mu).^2));
    
    % standardize the regressor variables
    ZR = 0.5+0.5*tanh((Z-mu)/Delta);
    % ZR = (Z-min(Z(isrf)))./(max(Z(isrf))-min(Z(isrf)));
    %
    Si2C = (aa*ZR + bb)./Z;   par.Si2C = Si2C;
    dSi2Cdaa = ZR./Z;         par.dSi2Cdaa = dSi2Cdaa;
    dSi2Cdbb = 1./Z;          par.dSi2Cdbb = dSi2Cdbb;
    %%%%%%%%%%%%%%%%%%%%%%% end %%%%%%%%%%%%%%%%%%%%%%%%

    TRdiv = par.TRdiv;
    PFdiv = buildPFD(par,'bSi');
    G     = uptake_Si(par);
    par.PFDs = PFdiv;
    % build Jacobian matrix
    % +++++++++++++++++++++++++++++++++++++++++
    % F = ([TRdiv*SIL + G*Si2C*SIL - kappa_si*DSI + kappa_gs*(SIL-SILbar);...
    % PFdiv*DSI - G*Si2C*SIL + kappa_si*DSI]);
    tic 
    Jac = [[TRdiv+kappa_gs*I+d0(G*Si2C), -d0(kappa_si)];...
           [                -d0(G*Si2C),  PFdiv+d0(kappa_si)]];
    
    RHS = [[SILbar*kappa_gs]; ...
           [sparse(nwet,1)]];
    
    FD = mfactor(Jac);
    Si = mfactor(FD,RHS);
    SIL = Si(1:nwet);
    DSI = Si(nwet+1:end);

    fprintf('Solving Si model...Elapsed time is %4.4e seconds \n',toc)
    %% +++++++++++++++++++++++++++++++++++++++++++
    if (par.optim == off)
        Six = [];
    else
        [~,Gx] = uptake_Si(par);
        % sigma
        if (par.opt_sigma == on)
            tmp =  [-d0(Gx(:,pindx.lsigma))*SIL.*Si2C;...
                    d0(Gx(:,pindx.lsigma))*SIL.*Si2C]   ;
            
            Six(:,pindx.lsigma) = mfactor(FD, tmp);
        end
        
        % bP
        if (par.opt_bP == on)
            tmp = [-d0(Gx(:,pindx.lbP))*SIL.*Si2C ;...
                   d0(Gx(:,pindx.lbP))*SIL.*Si2C]   ;
            
            Six(:,pindx.lbP) = mfactor(FD, tmp);
        end
        
        % bP_T
        if (par.opt_bP_T == on)
            tmp = [-d0(Gx(:,pindx.bP_T))*SIL.*Si2C ;...
                   d0(Gx(:,pindx.bP_T))*SIL.*Si2C];
            
            Six(:,pindx.bP_T) = mfactor(FD, tmp);
        end
        
        % kdP
        if (par.opt_kdP == on)
            tmp = [-d0(Gx(:,pindx.lkdP))*SIL.*Si2C;...
                   d0(Gx(:,pindx.lkdP))*SIL.*Si2C];
            
            Six(:,pindx.lkdP) = mfactor(FD, tmp);
        end
        
        % alpha
        if (par.opt_alpha == on)
            tmp = [-d0(Gx(:,pindx.lalpha))*SIL.*Si2C; ...
                   d0(Gx(:,pindx.lalpha))*SIL.*Si2C];
            
            Six(:,pindx.lalpha) = mfactor(FD, tmp);
        end
        
        % beta
        if (par.opt_beta == on)
            tmp = [-d0(Gx(:,pindx.lbeta))*SIL.*Si2C; ...
                   d0(Gx(:,pindx.lbeta))*SIL.*Si2C];
            
            Six(:,pindx.lbeta) = mfactor(FD, tmp);
        end
        
        % dsi
        if (par.opt_dsi == on)
            [~,Gout] = buildPFD(par,'bSi');
            PFD_dsi = Gout.PFD_d;
            par.PFD_dsi = PFD_dsi;
            tmp = dsi*[DSI*0;  -PFD_dsi*DSI];
            
            Six(:,pindx.ldsi) = mfactor(FD, tmp);
        end
        
        % at
        if (par.opt_at == on)
            [~,Gout] = buildPFD(par,'bSi');
            PFD_at = Gout.PFD_at;
            k_at = d0(exp(-bt./T));
            par.PFD_at = PFD_at;
            par.k_at = k_at;
            tmp = at*[k_at*DSI; ...
                      -(k_at+PFD_at)*DSI];
            
            Six(:,pindx.lat) = mfactor(FD, tmp);
        end
        
        % bt
        if (par.opt_bt == on)
            [~,Gout] = buildPFD(par,'bSi');
            PFD_bt  = Gout.PFD_bt;
            k_bt = -d0((at*exp(-bt./T))./T);
            par.PFD_bt = PFD_bt;
            par.k_bt   = k_bt;
            tmp = bt*[k_bt*DSI; ...
                      -(k_bt+PFD_bt)*DSI];
            
            Six(:,pindx.lbt) = mfactor(FD, tmp);
        end
        
        % aa
        if (par.opt_aa == on)
            tmp = [-d0(G*SIL)*dSi2Cdaa; ...
                   d0(G*SIL)*dSi2Cdaa];
            
            Six(:,pindx.aa) = mfactor(FD, tmp);
        end
        
        % bb
        if (par.opt_bb == on)
            tmp = bb*[-d0(G*SIL)*dSi2Cdbb; ...
                      d0(G*SIL)*dSi2Cdbb];
            
            Six(:,pindx.lbb) = mfactor(FD, tmp);
        end
    end
    SILx = Six(1:nwet,:);
    DSIx = Six(nwet+1:end,:);
    %% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    if (par.optim == off)
        Sixx = [];
    else 
        nx = par.npx + par.nsx;
        ncs = nchoosek(nx,2)+nx;
        Sixx = sparse(2*nwet,ncs);
        % Compute the hessian of the solution wrt the parameters
        DIP = par.DIP;
        DIPx  = par.Px(1:nwet,:);
        DIPxx = par.Pxx(1:nwet,:);
        [~,~,Gxx,Gp,par] = uptake_Si(par);
        dGpdalpha = par.dGpdalpha;
        dGpdbeta  = par.dGpdbeta;
        % sigma sigma
        kk = 1;
        if (par.opt_sigma)
            tmp = [-d0(Gxx(:,kk))*SIL.*Si2C; ...
                   d0(Gxx(:,kk))*SIL.*Si2C] + ...
                  2*[-d0(Gx(:,pindx.lsigma))*SILx(:,pindx.lsigma).*Si2C;...
                     d0(Gx(:,pindx.lsigma))*SILx(:,pindx.lsigma).*Si2C];
            
            Sixx(:,kk) = mfactor(FD, tmp);
            kk = kk + 1;
        end
        % sigma kdP
        if (par.opt_sigma & par.opt_kdP)
            tmp = [-d0(Gxx(:,kk))*SIL.*Si2C; ...
                   d0(Gxx(:,kk))*SIL.*Si2C] + ...
                  [-d0(Gx(:,pindx.lsigma))*SILx(:,pindx.lkdP).*Si2C;...
                   d0(Gx(:,pindx.lsigma))*SILx(:,pindx.lkdP).*Si2C] + ...
                  [-d0(Gx(:,pindx.lkdP))*SILx(:,pindx.lsigma).*Si2C;...
                   d0(Gx(:,pindx.lkdP))*SILx(:,pindx.lsigma).*Si2C];
            
            Sixx(:,kk) = mfactor(FD, tmp);
            kk = kk + 1;
        end
        % sigma bP_T
        if (par.opt_sigma & par.opt_bP_T)
            tmp = [-d0(Gxx(:,kk))*SIL.*Si2C; ...
                   d0(Gxx(:,kk))*SIL.*Si2C] + ...
                  [-d0(Gx(:,pindx.lsigma))*SILx(:,pindx.bP_T).*Si2C;...
                   d0(Gx(:,pindx.lsigma))*SILx(:,pindx.bP_T).*Si2C] + ...
                  [-d0(Gx(:,pindx.bP_T))*SILx(:,pindx.lsigma).*Si2C;...
                   d0(Gx(:,pindx.bP_T))*SILx(:,pindx.lsigma).*Si2C];
            
            Sixx(:,kk) = mfactor(FD, tmp);
            kk = kk + 1;
        end
        % sigma bP
        if (par.opt_sigma == on & par.opt_bP == on)
            tmp = [-d0(Gxx(:,kk))*SIL.*Si2C; ...
                   d0(Gxx(:,kk))*SIL.*Si2C] + ...
                  [-d0(Gx(:,pindx.lsigma))*SILx(:,pindx.lbP).*Si2C;...
                   d0(Gx(:,pindx.lsigma))*SILx(:,pindx.lbP).*Si2C]+...
                  [-d0(Gx(:,pindx.lbP))*SILx(:,pindx.lsigma).*Si2C;...
                   d0(Gx(:,pindx.lbP))*SILx(:,pindx.lsigma).*Si2C];
            
            Sixx(:,kk) = mfactor(FD, tmp);
            kk = kk + 1;
        end
        % sigma alpha
        if (par.opt_sigma == on & par.opt_alpha == on)
            tmp = [-d0(Gxx(:,kk))*SIL.*Si2C; ...
                   d0(Gxx(:,kk))*SIL.*Si2C] + ...
                  [-d0(Gx(:,pindx.lsigma))*SILx(:,pindx.lalpha).*Si2C;...
                   d0(Gx(:,pindx.lsigma))*SILx(:,pindx.lalpha).*Si2C]+...
                  [-d0(Gx(:,pindx.lalpha))*SILx(:,pindx.lsigma).*Si2C;...
                   d0(Gx(:,pindx.lalpha))*SILx(:,pindx.lsigma).*Si2C];
            
            Sixx(:,kk) = mfactor(FD, tmp);
            kk = kk + 1;
        end
        % sigma beta
        if (par.opt_sigma == on & par.opt_beta == on)
            tmp = [-d0(Gxx(:,kk))*SIL.*Si2C; ...
                   d0(Gxx(:,kk))*SIL.*Si2C] + ...
                  [-d0(Gx(:,pindx.lsigma))*SILx(:,pindx.lbeta).*Si2C;...
                   d0(Gx(:,pindx.lsigma))*SILx(:,pindx.lbeta).*Si2C]+...
                  [-d0(Gx(:,pindx.lbeta))*SILx(:,pindx.lsigma).*Si2C;...
                   d0(Gx(:,pindx.lbeta))*SILx(:,pindx.lsigma).*Si2C];
            
            Sixx(:,kk) = mfactor(FD, tmp);
            kk = kk + 1;
        end
        % kdP kdP
        if (par.opt_kdP == on)
            tmp = [-d0(Gxx(:,kk))*SIL.*Si2C; ...
                   d0(Gxx(:,kk))*SIL.*Si2C] + ...
                  2*[-d0(Gx(:,pindx.lkdP))*SILx(:,pindx.lkdP).*Si2C;...
                     d0(Gx(:,pindx.lkdP))*SILx(:,pindx.lkdP).*Si2C];
            
            Sixx(:,kk) = mfactor(FD, tmp);
            kk = kk + 1;
        end
        % kdP bP_T
        if (par.opt_kdP == on & par.opt_bP_T == on)
            tmp = [-d0(Gxx(:,kk))*SIL.*Si2C; ...
                   d0(Gxx(:,kk))*SIL.*Si2C] + ...
                  [-d0(Gx(:,pindx.lkdP))*SILx(:,pindx.bP_T).*Si2C;...
                   d0(Gx(:,pindx.lkdP))*SILx(:,pindx.bP_T).*Si2C]+...
                  [-d0(Gx(:,pindx.bP_T))*SILx(:,pindx.lkdP).*Si2C;...
                   d0(Gx(:,pindx.bP_T))*SILx(:,pindx.lkdP).*Si2C];
            
            Sixx(:,kk) = mfactor(FD, tmp);
            kk = kk + 1;
        end
        % kdP bP
        if (par.opt_kdP == on & par.opt_bP == on)
            tmp = [-d0(Gxx(:,kk))*SIL.*Si2C; ...
                   d0(Gxx(:,kk))*SIL.*Si2C] + ...
                  [-d0(Gx(:,pindx.lkdP))*SILx(:,pindx.lbP).*Si2C;...
                   d0(Gx(:,pindx.lkdP))*SILx(:,pindx.lbP).*Si2C]+...
                  [-d0(Gx(:,pindx.lbP))*SILx(:,pindx.lkdP).*Si2C;...
                   d0(Gx(:,pindx.lbP))*SILx(:,pindx.lkdP).*Si2C];
            
            Sixx(:,kk) = mfactor(FD, tmp);
            kk = kk + 1;
        end
        % kdP alpha
        if (par.opt_kdP == on & par.opt_alpha == on)
            tmp = [-d0(Gxx(:,kk))*SIL.*Si2C; ...
                   d0(Gxx(:,kk))*SIL.*Si2C] + ...
                  [-d0(Gx(:,pindx.lkdP))*SILx(:,pindx.lalpha).*Si2C;...
                   d0(Gx(:,pindx.lkdP))*SILx(:,pindx.lalpha).*Si2C]+...
                  [-d0(Gx(:,pindx.lalpha))*SILx(:,pindx.lkdP).*Si2C;...
                   d0(Gx(:,pindx.lalpha))*SILx(:,pindx.lkdP).*Si2C];
            
            Sixx(:,kk) = mfactor(FD, tmp);
            kk = kk + 1;
        end
        % kdP beta
        if (par.opt_kdP == on & par.opt_beta == on)
            tmp = [-d0(Gxx(:,kk))*SIL.*Si2C; ...
                   d0(Gxx(:,kk))*SIL.*Si2C] + ...
                  [-d0(Gx(:,pindx.lkdP))*SILx(:,pindx.lbeta).*Si2C;...
                   d0(Gx(:,pindx.lkdP))*SILx(:,pindx.lbeta).*Si2C]+...
                  [-d0(Gx(:,pindx.lbeta))*SILx(:,pindx.lkdP).*Si2C;...
                   d0(Gx(:,pindx.lbeta))*SILx(:,pindx.lkdP).*Si2C];
            
            Sixx(:,kk) = mfactor(FD, tmp);
            kk = kk + 1;
        end
        % bP_T bP_T
        if (par.opt_bP_T == on)
            tmp = [-d0(Gxx(:,kk))*SIL.*Si2C; ...
                   d0(Gxx(:,kk))*SIL.*Si2C] + ...
                  2*[-d0(Gx(:,pindx.bP_T))*SILx(:,pindx.bP_T).*Si2C;...
                     d0(Gx(:,pindx.bP_T))*SILx(:,pindx.bP_T).*Si2C];
            
            Sixx(:,kk) = mfactor(FD, tmp);
            kk = kk + 1;
        end
        % bP_T bP
        if (par.opt_bP_T == on & par.opt_bP == on)
            tmp = [-d0(Gxx(:,kk))*SIL.*Si2C; ...
                   d0(Gxx(:,kk))*SIL.*Si2C] + ...
                  [-d0(Gx(:,pindx.bP_T))*SILx(:,pindx.lbP).*Si2C;...
                   d0(Gx(:,pindx.bP_T))*SILx(:,pindx.lbP).*Si2C]+...
                  [-d0(Gx(:,pindx.lbP))*SILx(:,pindx.bP_T).*Si2C;...
                   d0(Gx(:,pindx.lbP))*SILx(:,pindx.bP_T).*Si2C];
            
            Sixx(:,kk) = mfactor(FD, tmp);
            kk = kk + 1;
        end
        % bP_T alpha
        if (par.opt_bP_T == on & par.opt_alpha == on)
            tmp = [-d0(Gxx(:,kk))*SIL.*Si2C; ...
                   d0(Gxx(:,kk))*SIL.*Si2C] + ...
                  [-d0(Gx(:,pindx.bP_T))*SILx(:,pindx.lalpha).*Si2C;...
                   d0(Gx(:,pindx.bP_T))*SILx(:,pindx.lalpha).*Si2C]+...
                  [-d0(Gx(:,pindx.lalpha))*SILx(:,pindx.bP_T).*Si2C;...
                   d0(Gx(:,pindx.lalpha))*SILx(:,pindx.bP_T).*Si2C];
            
            Sixx(:,kk) = mfactor(FD, tmp);
            kk = kk + 1;
        end
        % bP_T beta
        if (par.opt_bP_T == on & par.opt_beta == on)
            tmp = [-d0(Gxx(:,kk))*SIL.*Si2C; ...
                   d0(Gxx(:,kk))*SIL.*Si2C] + ...
                  [-d0(Gx(:,pindx.bP_T))*SILx(:,pindx.lbeta).*Si2C;...
                   d0(Gx(:,pindx.bP_T))*SILx(:,pindx.lbeta).*Si2C]+...
                  [-d0(Gx(:,pindx.lbeta))*SILx(:,pindx.bP_T).*Si2C;...
                   d0(Gx(:,pindx.lbeta))*SILx(:,pindx.bP_T).*Si2C];
            
            Sixx(:,kk) = mfactor(FD, tmp);
            kk = kk + 1;
        end
        % bP bP
        if (par.opt_bP == on)
            tmp = [-d0(Gxx(:,kk))*SIL.*Si2C; ...
                   d0(Gxx(:,kk))*SIL.*Si2C] + ...
                  2*[-d0(Gx(:,pindx.lbP))*SILx(:,pindx.lbP).*Si2C;...
                     d0(Gx(:,pindx.lbP))*SILx(:,pindx.lbP).*Si2C];
            
            Sixx(:,kk) = mfactor(FD, tmp);
            kk = kk + 1;
        end
        % bP alpha
        if (par.opt_bP == on & par.opt_alpha == on)
            tmp = [-d0(Gxx(:,kk))*SIL.*Si2C; ...
                   d0(Gxx(:,kk))*SIL.*Si2C] + ...
                  [-d0(Gx(:,pindx.lbP))*SILx(:,pindx.lalpha).*Si2C; ...
                   d0(Gx(:,pindx.lbP))*SILx(:,pindx.lalpha).*Si2C]+ ...
                  [-d0(Gx(:,pindx.lalpha))*SILx(:,pindx.lbP).*Si2C; ...
                   d0(Gx(:,pindx.lalpha))*SILx(:,pindx.lbP).*Si2C];
            
            Sixx(:,kk) = mfactor(FD, tmp);
            kk = kk + 1;
        end
        % bP beta
        if (par.opt_bP == on & par.opt_beta == on)
            tmp = [-d0(Gxx(:,kk))*SIL.*Si2C; ...
                   d0(Gxx(:,kk))*SIL.*Si2C] + ...
                  [-d0(Gx(:,pindx.lbP))*SILx(:,pindx.lbeta).*Si2C; ...
                   d0(Gx(:,pindx.lbP))*SILx(:,pindx.lbeta).*Si2C]+ ...
                  [-d0(Gx(:,pindx.lbeta))*SILx(:,pindx.lbP).*Si2C; ...
                   d0(Gx(:,pindx.lbeta))*SILx(:,pindx.lbP).*Si2C];
            
            Sixx(:,kk) = mfactor(FD, tmp);
            kk = kk + 1;
        end
        % alpha alpha
        if (par.opt_alpha == on)
            tmp = [-d0(Gxx(:,kk))*SIL.*Si2C; ...
                   d0(Gxx(:,kk))*SIL.*Si2C] + ...
                  2*[-d0(Gx(:,pindx.lalpha))*SILx(:,pindx.lalpha).*Si2C; ...
                     d0(Gx(:,pindx.lalpha))*SILx(:,pindx.lalpha).*Si2C];
            
            Sixx(:,kk) = mfactor(FD, tmp);
            kk = kk + 1;
        end
        % alpha beta
        if (par.opt_alpha == on & par.opt_beta == on)
            tmp = [-d0(Gxx(:,kk))*SIL.*Si2C; ...
                   d0(Gxx(:,kk))*SIL.*Si2C] + ...
                  [-d0(Gx(:,pindx.lalpha))*SILx(:,pindx.lbeta).*Si2C; ...
                   d0(Gx(:,pindx.lalpha))*SILx(:,pindx.lbeta).*Si2C]+ ...
                  [-d0(Gx(:,pindx.lbeta))*SILx(:,pindx.lalpha).*Si2C; ...
                   d0(Gx(:,pindx.lbeta))*SILx(:,pindx.lalpha).*Si2C];
            
            Sixx(:,kk) = mfactor(FD, tmp);
            kk = kk + 1;
        end
        % beta beta
        if (par.opt_beta == on)
            tmp = [-d0(Gxx(:,kk))*SIL.*Si2C; ...
                   d0(Gxx(:,kk))*SIL.*Si2C] + ...
                  2*[-d0(Gx(:,pindx.lbeta))*SILx(:,pindx.lbeta).*Si2C; ...
                     d0(Gx(:,pindx.lbeta))*SILx(:,pindx.lbeta).*Si2C];
            
            Sixx(:,kk) = mfactor(FD, tmp);
            kk = kk + 1;
        end
        % --------------------------------------------------------
        % Simodel parameters
        % sigma dsi
        if (par.opt_sigma == on & par.opt_dsi == on)
            tmp = dsi*[0*DSI; ...
                       -PFD_dsi*DSIx(:,pindx.lsigma)] + ...
                  [-d0(Gx(:,pindx.lsigma))*SILx(:,pindx.ldsi).*Si2C; ...
                   d0(Gx(:,pindx.lsigma))*SILx(:,pindx.ldsi).*Si2C];
            
            Sixx(:,kk) = mfactor(FD, tmp);
            kk = kk + 1;
        end
        % sigma at
        if (par.opt_sigma == on & par.opt_at == on)
            tmp = at*[k_at*DSIx(:,pindx.lsigma); ...
                      -(k_at + PFD_at)*DSIx(:,pindx.lsigma)]+ ...
                  [-d0(Gx(:,pindx.lsigma))*SILx(:,pindx.lat).*Si2C; ...
                   d0(Gx(:,pindx.lsigma))*SILx(:,pindx.lat).*Si2C];
            
            Sixx(:,kk) = mfactor(FD, tmp);
            kk = kk + 1;
        end
        % sigma bt
        if (par.opt_sigma == on & par.opt_bt == on)
            tmp = bt*[k_bt*DSIx(:,pindx.lsigma); ...
                      -(k_bt+PFD_bt)*DSIx(:,pindx.lsigma)]+ ...
                  [-d0(Gx(:,pindx.lsigma))*SILx(:,pindx.lbt).*Si2C; ...
                   d0(Gx(:,pindx.lsigma))*SILx(:,pindx.lbt).*Si2C];
            
            Sixx(:,kk) = mfactor(FD, tmp);
            kk = kk + 1;
        end
        % sigma aa
        if (par.opt_sigma == on & par.opt_aa == on)
            tmp = [-d0(Gxx(:,kk))*SIL; ...
                   d0(Gxx(:,kk))*SIL] + ...
                  [-d0(G*SILx(:,pindx.lsigma))*dSi2Cdaa; ...
                   d0(G*SILx(:,pindx.lsigma))*dSi2Cdaa]+ ...
                  [-d0(Gx(:,pindx.lsigma))*SILx(:,pindx.aa).*Si2C; ...
                   d0(Gx(:,pindx.lsigma))*SILx(:,pindx.aa).*Si2C];
            
            Sixx(:,kk) = mfactor(FD, tmp);
            kk = kk + 1;
        end
        % sigma bb
        if (par.opt_sigma == on & par.opt_bb == on)
            tmp = [-d0(Gxx(:,kk))*SIL; ...
                   d0(Gxx(:,kk))*SIL] + ...
                  bb*[-d0(G*SILx(:,pindx.lsigma))*dSi2Cdbb; ...
                      d0(G*SILx(:,pindx.lsigma))*dSi2Cdbb]+ ...
                  [-d0(Gx(:,pindx.lsigma))*SILx(:,pindx.lbb).*Si2C; ...
                   d0(Gx(:,pindx.lsigma))*SILx(:,pindx.lbb).*Si2C];
            
            Sixx(:,kk) = mfactor(FD, tmp);
            kk = kk + 1;
        end
        % kdP dsi
        if (par.opt_kdP == on & par.opt_dsi == on)
            tmp = dsi*[0*DSI; ...
                       -(PFD_dsi)*DSIx(:,pindx.lkdP)]+ ...
                  [-d0(Gx(:,pindx.lkdP))*SILx(:,pindx.ldsi).*Si2C; ...
                   d0(Gx(:,pindx.lkdP))*SILx(:,pindx.ldsi).*Si2C];

            Sixx(:,kk) = mfactor(FD, tmp);
            kk = kk + 1;
        end
        % kdP at
        if (par.opt_kdP == on & par.opt_at == on)
            tmp = at*[k_at*DSIx(:,pindx.lkdP); ...
                      -(k_at+PFD_at)*DSIx(:,pindx.lkdP)]+ ...
                  [-d0(Gx(:,pindx.lkdP))*SILx(:,pindx.lat).*Si2C; ...
                   d0(Gx(:,pindx.lkdP))*SILx(:,pindx.lat).*Si2C];

            Sixx(:,kk) = mfactor(FD, tmp);
            kk = kk + 1;
        end
        % kdP bt
        if (par.opt_kdP == on & par.opt_bt == on)
            tmp = bt*[k_bt*DSIx(:,pindx.lkdP); ...
                      -(k_bt+PFD_bt)*DSIx(:,pindx.lkdP)]+ ...
                  [-d0(Gx(:,pindx.lkdP))*SILx(:,pindx.lbt).*Si2C; ...
                   d0(Gx(:,pindx.lkdP))*SILx(:,pindx.lbt).*Si2C];

            Sixx(:,kk) = mfactor(FD, tmp); 
            kk = kk + 1;
        end
        % kdP aa
        if (par.opt_kdP == on & par.opt_aa == on)
            tmp =  [-d0(Gxx(:,kk))*SIL; ...
                    d0(Gxx(:,kk))*SIL] + ...
                   [-d0(G*SILx(:,pindx.lkdP))*dSi2Cdaa; ...
                    d0(G*SILx(:,pindx.lkdP))*dSi2Cdaa]+ ...
                   [-d0(Gx(:,pindx.lkdP))*SILx(:,pindx.aa).*Si2C; ...
                    d0(Gx(:,pindx.lkdP))*SILx(:,pindx.aa).*Si2C];

            Sixx(:,kk) = mfactor(FD, tmp);
            kk = kk + 1;
        end
        % kdP bb
        if (par.opt_kdP == on & par.opt_bb == on)
            tmp = [-d0(Gxx(:,kk))*SIL; ...
                   d0(Gxx(:,kk))*SIL] + ...
                  bb*[-d0(G*SILx(:,pindx.lkdP))*dSi2Cdbb; ...
                      d0(G*SILx(:,pindx.lkdP))*dSi2Cdbb]+ ...
                  [-d0(Gx(:,pindx.lkdP))*SILx(:,pindx.lbb).*Si2C; ...
                   d0(Gx(:,pindx.lkdP))*SILx(:,pindx.lbb).*Si2C];

            Sixx(:,kk) = mfactor(FD, tmp);
            kk = kk + 1;
        end
        % slope dsi
        if (par.opt_bP_T == on & par.opt_dsi == on)
            tmp = dsi*[0*DSI; ...
                       -PFD_dsi*DSIx(:,pindx.bP_T)]+ ...
                  [-d0(Gx(:,pindx.bP_T))*SILx(:,pindx.ldsi).*Si2C; ...
                   d0(Gx(:,pindx.bP_T))*SILx(:,pindx.ldsi).*Si2C];

            Sixx(:,kk) = mfactor(FD, tmp);
            kk = kk + 1;
        end
        % slope at
        if (par.opt_bP_T == on & par.opt_at == on)
            tmp = at*[k_at*DSIx(:,pindx.bP_T); ...
                      -(k_at+PFD_at)*DSIx(:,pindx.bP_T)]+ ...
                  [-d0(Gx(:,pindx.bP_T))*SILx(:,pindx.lat).*Si2C; ...
                   d0(Gx(:,pindx.bP_T))*SILx(:,pindx.lat).*Si2C];

            Sixx(:,kk) = mfactor(FD, tmp);
            kk = kk + 1;
        end
        % slope bt
        if (par.opt_bP_T == on & par.opt_bt == on)
            tmp = bt*[k_bt*DSIx(:,pindx.bP_T); ...
                      -(k_bt+PFD_bt)*DSIx(:,pindx.bP_T)]+ ...
                  [-d0(Gx(:,pindx.bP_T))*SILx(:,pindx.lbt).*Si2C; ...
                   d0(Gx(:,pindx.bP_T))*SILx(:,pindx.lbt).*Si2C];

            Sixx(:,kk) = mfactor(FD, tmp);
            kk = kk + 1;
        end
        % bP_T aa
        if (par.opt_bP_T == on & par.opt_aa == on)
            tmp = [-d0(Gxx(:,kk))*SIL; ...
                   d0(Gxx(:,kk))*SIL] + ...
                  [-d0(G*SILx(:,pindx.bP_T))*dSi2Cdaa; ...
                   d0(G*SILx(:,pindx.bP_T))*dSi2Cdaa] + ...
                  [-d0(Gx(:,pindx.bP_T))*SILx(:,pindx.aa).*Si2C;  ...
                   d0(Gx(:,pindx.bP_T))*SILx(:,pindx.aa).*Si2C];

            Sixx(:,kk) = mfactor(FD, tmp);
            kk = kk + 1;
        end
        % bP_T bb
        if (par.opt_bP_T == on & par.opt_bb == on)
            tmp = [-d0(Gxx(:,kk))*SIL; ...
                   d0(Gxx(:,kk))*SIL] + ...
                  bb*[-d0(G*SILx(:,pindx.bP_T))*dSi2Cdbb; ...
                      d0(G*SILx(:,pindx.bP_T))*dSi2Cdbb] + ...
                  [-d0(Gx(:,pindx.bP_T))*SILx(:,pindx.lbb).*Si2C; ...
                   d0(Gx(:,pindx.bP_T))*SILx(:,pindx.lbb).*Si2C];

            Sixx(:,kk) = mfactor(FD, tmp);
            kk = kk + 1;
        end
        % bP dsi
        if (par.opt_bP == on & par.opt_dsi == on)
            tmp = dsi*[0*DSI; ...
                       -PFD_dsi*DSIx(:, pindx.lbP)] + ...
                  [-d0(Gx(:,pindx.lbP))*SILx(:,pindx.ldsi).*Si2C; ...
                   d0(Gx(:,pindx.lbP))*SILx(:,pindx.ldsi).*Si2C];

            Sixx(:,kk) = mfactor(FD, tmp);
            kk = kk + 1;
        end
        % bP at
        if (par.opt_bP == on & par.opt_at == on)
            tmp = at*[k_at*DSIx(:, pindx.lbP); ...
                      -(k_at+PFD_at)*DSIx(:, pindx.lbP)]+ ...
                  [-d0(Gx(:,pindx.lbP))*SILx(:,pindx.lat).*Si2C; ...
                   d0(Gx(:,pindx.lbP))*SILx(:,pindx.lat).*Si2C];

            Sixx(:,kk) = mfactor(FD, tmp);
            kk = kk + 1;
        end
        % bP bt
        if (par.opt_bP == on & par.opt_bt == on)
            tmp = bt*[k_bt*DSIx(:, pindx.lbP); ...
                      -(k_bt+PFD_bt)*DSIx(:,pindx.lbP)]+ ...
                  [-d0(Gp*DIPx(:,pindx.lbP))*SILx(:,pindx.lbt).*Si2C; ...
                   d0(Gp*DIPx(:,pindx.lbP))*SILx(:,pindx.lbt).*Si2C];

            Sixx(:,kk) = mfactor(FD, tmp);
            kk = kk + 1;
        end
        % bP aa
        if (par.opt_bP == on & par.opt_aa == on)
            tmp = [-d0(Gxx(:,kk))*SIL; ...
                   d0(Gxx(:,kk))*SIL] + ...
                  [-d0(G*SILx(:,pindx.lbP))*dSi2Cdaa; ...
                   d0(G*SILx(:,pindx.lbP))*dSi2Cdaa] + ...
                  [-d0(Gx(:,pindx.lbP))*SILx(:,pindx.aa).*Si2C; ...
                   d0(Gx(:,pindx.lbP))*SILx(:,pindx.aa).*Si2C];

            Sixx(:,kk) = mfactor(FD, tmp);
            kk = kk + 1;
        end
        % bP bb
        if (par.opt_bP == on & par.opt_bb == on)
            tmp = [-d0(Gxx(:,kk))*SIL; ...
                   d0(Gxx(:,kk))*SIL] + ...
                  bb*[-d0(G*SILx(:,pindx.lbP))*dSi2Cdbb; ...
                      d0(G*SILx(:,pindx.lbP))*dSi2Cdbb]+ ...
                  [-d0(Gx(:,pindx.lbP))*SILx(:,pindx.lbb).*Si2C; ...
                   d0(Gx(:,pindx.lbP))*SILx(:,pindx.lbb).*Si2C];

            Sixx(:,kk) = mfactor(FD, tmp);
            kk = kk + 1;
        end
        % alpha dsi
        if (par.opt_alpha == on & par.opt_dsi == on)
            tmp = dsi*[0*DSI; ...
                       -PFD_dsi*DSIx(:,pindx.lalpha)] + ...
                  [-d0(Gx(:,pindx.lalpha))*SILx(:,pindx.ldsi).*Si2C; ...
                   d0(Gx(:,pindx.lalpha))*SILx(:,pindx.ldsi).*Si2C];

            Sixx(:,kk) = mfactor(FD, tmp);
            kk = kk + 1;
        end
        % alpha at
        if (par.opt_alpha == on & par.opt_at == on)
            tmp = at*[k_at*DSIx(:,pindx.lalpha); ...
                      -(k_at+PFD_at)*DSIx(:,pindx.lalpha)] + ...
                  [-d0(Gx(:,pindx.lalpha))*SILx(:,pindx.lat).*Si2C; ...
                   d0(Gx(:,pindx.lalpha))*SILx(:,pindx.lat).*Si2C];

            Sixx(:,kk) = mfactor(FD, tmp);
            kk = kk + 1;
        end
        % alpha bt
        if (par.opt_alpha == on & par.opt_bt == on)
            tmp = bt*[k_bt*DSIx(:,pindx.lalpha); ...
                      -(k_bt + PFD_bt)*DSIx(:,pindx.lalpha)] + ...
                  [-d0(Gx(:,pindx.lalpha))*SILx(:,pindx.lbt).*Si2C; ...
                   d0(Gx(:,pindx.lalpha))*SILx(:,pindx.lbt).*Si2C];

            Sixx(:,kk) = mfactor(FD, tmp);
            kk = kk + 1;
        end
        % alpha aa
        if (par.opt_alpha == on & par.opt_aa == on)
            tmp = [-d0(Gxx(:,kk))*SIL; ...
                   d0(Gxx(:,kk))*SIL] + ...
                  [-d0(G*SILx(:,pindx.lalpha))*dSi2Cdaa; ...
                   d0(G*SILx(:,pindx.lalpha))*dSi2Cdaa] + ...
                  [-d0(Gx(:,pindx.lalpha))*SILx(:,pindx.aa).*Si2C; ...
                   d0(Gx(:,pindx.lalpha))*SILx(:,pindx.aa).*Si2C];

            Sixx(:,kk) = mfactor(FD, tmp);
            kk = kk + 1;
        end
        % alpha bb
        if (par.opt_alpha == on & par.opt_bb == on)
            tmp = [-d0(Gxx(:,kk))*SIL; ...
                   d0(Gxx(:,kk))*SIL] + ...
                  bb*[-d0(G*SILx(:,pindx.lalpha))*dSi2Cdbb; ...
                      d0(G*SILx(:,pindx.lalpha))*dSi2Cdbb] + ...
                  [-d0(Gx(:,pindx.lalpha))*SILx(:,pindx.lbb).*Si2C; ...
                   d0(Gx(:,pindx.lalpha))*SILx(:,pindx.lbb).*Si2C];

            Sixx(:,kk) = mfactor(FD, tmp);
            kk = kk + 1;
        end
        % beta dsi
        if (par.opt_beta == on & par.opt_dsi == on)
            tmp = dsi*[0*DSI; ...
                       -PFD_dsi*DSIx(:,pindx.lbeta)] + ...
                  [-d0(Gx(:,pindx.lbeta))*SILx(:,pindx.ldsi).*Si2C; ...
                   d0(Gx(:,pindx.lbeta))*SILx(:,pindx.ldsi).*Si2C];
            Sixx(:,kk) = mfactor(FD, tmp);
            kk = kk + 1;
        end
        % beta at
        if (par.opt_beta == on & par.opt_at == on)
            tmp = at*[k_at*DSIx(:,pindx.lbeta); ...
                      -(k_at+PFD_at)*DSIx(:,pindx.lbeta)] + ...
                  [-d0(Gx(:,pindx.lbeta))*SILx(:,pindx.lat).*Si2C; ...
                   d0(Gx(:,pindx.lbeta))*SILx(:,pindx.lat).*Si2C];

            Sixx(:,kk) = mfactor(FD, tmp);
            kk = kk + 1;
        end
        % beta bt
        if (par.opt_beta == on & par.opt_bt == on)
            tmp = bt*[k_bt*DSIx(:,pindx.lbeta); ...
                      -(k_bt+PFD_bt)*DSIx(:,pindx.lbeta)] + ...
                  [-d0(Gx(:,pindx.lbeta))*SILx(:,pindx.lbt).*Si2C; ...
                   d0(Gx(:,pindx.lbeta))*SILx(:,pindx.lbt).*Si2C];

            Sixx(:,kk) = mfactor(FD, tmp);
            kk = kk + 1;
        end
        % beta aa
        if (par.opt_beta == on & par.opt_aa == on)
            tmp = [-d0(Gxx(:,kk))*SIL; ...
                   d0(Gxx(:,kk))*SIL] + ...
                  [-d0(G*SILx(:,pindx.lbeta))*dSi2Cdaa; ...
                   d0(G*SILx(:,pindx.lbeta))*dSi2Cdaa] + ...
                  [-d0(Gx(:,pindx.lbeta))*SILx(:,pindx.aa).*Si2C; ...
                   d0(Gx(:,pindx.lbeta))*SILx(:,pindx.aa).*Si2C];

            Sixx(:,kk) = mfactor(FD, tmp);
            kk = kk + 1;
        end
        % beta bb
        if (par.opt_beta == on & par.opt_bb == on)
            tmp = [-d0(Gxx(:,kk))*SIL; ...
                   d0(Gxx(:,kk))*SIL] + ...
                  bb*[-d0(G*SILx(:,pindx.lbeta))*dSi2Cdbb; ...
                      d0(G*SILx(:,pindx.lbeta))*dSi2Cdbb] + ...
                  [-d0(Gx(:,pindx.lbeta))*SILx(:,pindx.lbb).*Si2C; ...
                   d0(Gx(:,pindx.lbeta))*SILx(:,pindx.lbb).*Si2C];

            Sixx(:,kk) = mfactor(FD, tmp);
            kk = kk + 1;
        end
        % dsi dsi
        if (par.opt_dsi == on)
            [~,~,Hout] = buildPFD(par,'bSi');
            PFD_dsi_dsi = Hout.PFD_d_d;
            par.PFD_dsi_dsi = PFD_dsi_dsi;
            tmp = dsi*[0*DSI; ...
                       -PFD_dsi*DSI] + ...
                  dsi*dsi*[0*DSI; ...
                           -PFD_dsi_dsi*DSI] + ...
                  2*dsi*[0*DSI; ...
                         -PFD_dsi*DSIx(:,pindx.ldsi)];

            Sixx(:,kk) = mfactor(FD, tmp);
            kk = kk + 1;
        end
        % dsi at
        if (par.opt_dsi & par.opt_at)
            [~,~,Hout] = buildPFD(par,'bSi');
            PFD_at_d = Hout.PFD_at_d;
            par.PFD_at_d = PFD_at_d;
            tmp = dsi*at*[0*DSI; ...
                          -PFD_at_d*DSI] + ...
                  dsi*[0*DSI; ...
                       -PFD_dsi*DSIx(:,pindx.lat)] + ...
                  at*[k_at*DSIx(:,pindx.ldsi); ...
                      -(k_at+PFD_at)*DSIx(:,pindx.ldsi)];
            
            Sixx(:,kk) = mfactor(FD, tmp);
            kk = kk + 1;
        end
        % dsi bt
        if (par.opt_dsi & par.opt_bt)
            [~,~,Hout] = buildPFD(par,'bSi');
            PFD_bt_dsi = Hout.PFD_bt_d;
            par.PFD_bt_dsi = PFD_bt_dsi;
            tmp = dsi*bt*[0*DSI; ...
                          -PFD_bt_dsi*DSI] + ...
                  dsi*[0*DSI; ...
                       -PFD_dsi*DSIx(:,pindx.lbt)] + ...
                  bt*[k_bt*DSIx(:,pindx.ldsi); ...
                      -(k_bt+PFD_bt)*DSIx(:,pindx.ldsi)];
            
            Sixx(:,kk) = mfactor(FD, tmp);
            kk = kk + 1;
        end
        % dsi aa
        if (par.opt_dsi & par.opt_aa)
            tmp = dsi*[0*DSI; ...
                       -PFD_dsi*DSIx(:,pindx.aa)] + ...
                  [-d0(G*SILx(:,pindx.ldsi))*dSi2Cdaa; ...
                   d0(G*SILx(:,pindx.ldsi))*dSi2Cdaa] + ...
                  [-d0(Gx(:,pindx.ldsi))*SILx(:,pindx.aa).*Si2C; ...
                   d0(Gx(:,pindx.ldsi))*SILx(:,pindx.aa).*Si2C];
            
            Sixx(:,kk) = mfactor(FD, tmp);
            kk = kk + 1;
        end
        % dsi bb
        if (par.opt_dsi & par.opt_bb)
            tmp = dsi*[0*DSI; ...
                       -PFD_dsi*DSIx(:,pindx.lbb)] + ...
                  bb*[-d0(G*SILx(:,pindx.ldsi))*dSi2Cdbb; ...
                      d0(G*SILx(:,pindx.ldsi))*dSi2Cdbb] + ...
                  [-d0(Gx(:,pindx.ldsi))*SILx(:,pindx.lbb).*Si2C; ...
                   d0(Gx(:,pindx.ldsi))*SILx(:,pindx.lbb).*Si2C];

            Sixx(:,kk) = mfactor(FD, tmp);
            kk = kk + 1;
        end
        % at at
        if (par.opt_at == on)
            [~,~,Hout] = buildPFD(par,'bSi');
            PFD_at_at  = Hout.PFD_at_at;
            k_at_at = 0;
            par.PFD_at_at = PFD_at_at;
            par.k_at_at = k_at_at;
            tmp = at*[k_at*DSI; ...
                      -(k_at+PFD_at)*DSI] + ...
                  at*at*[0*DSI; ...
                         -PFD_at_at*DSI] + ...
                  2*at*[k_at*DSIx(:,pindx.lat); ...
                        -(k_at+PFD_at)*DSIx(:,pindx.lat)];

            Sixx(:,kk) = mfactor(FD, tmp);
            kk = kk + 1;
        end
        % at bt
        if (par.opt_at == on & par.opt_bt == on)
            [~,~,Hout]  = buildPFD(par,'bSi');
            PFD_at_bt = Hout.PFD_at_bt;
            k_at_bt = -d0(exp(-bt./T).*(1./T));
            par.k_at_bt = k_at_bt;
            par.PFD_at_bt = PFD_at_bt;
            tmp = at*bt*[k_at_bt; ...
                         -(k_at_bt+PFD_at_bt)]*DSI + ...
                  at*[k_at*DSIx(:,pindx.lbt); ...
                      -(k_at+PFD_at)*DSIx(:,pindx.lbt)]+ ...
                  bt*[k_bt*DSIx(:,pindx.lat); ...
                      -(k_bt+PFD_bt)*DSIx(:,pindx.lat)];
            
            Sixx(:,kk) = mfactor(FD, tmp);
            kk = kk + 1;
        end
        % at aa
        if (par.opt_at == on & par.opt_aa == on)
            tmp = at*[k_at*DSIx(:,pindx.aa); ...
                      -(k_at+PFD_at)*DSIx(:,pindx.aa)] + ...
                  [-d0(G*SILx(:,pindx.lat))*dSi2Cdaa; ... 
                   d0(G*SILx(:,pindx.lat))*dSi2Cdaa] + ...
                  [-d0(Gx(:,pindx.lat))*SILx(:,pindx.aa).*Si2C; ...
                   d0(Gx(:,pindx.lat))*SILx(:,pindx.aa).*Si2C];
            
            Sixx(:,kk) = mfactor(FD, tmp);
            kk = kk + 1;
        end
        % at bb
        if (par.opt_at == on & par.opt_bb == on)
            tmp = at*[k_at*DSIx(:,pindx.lbb); ...
                      -(k_at+PFD_at)*DSIx(:,pindx.lbb)] + ...
                  bb*[-d0(G*SILx(:,pindx.lat))*dSi2Cdbb; ... 
                      d0(G*SILx(:,pindx.lat))*dSi2Cdbb] + ...
                  [-d0(Gx(:,pindx.lat))*SILx(:,pindx.lbb).*Si2C; ...
                   d0(Gx(:,pindx.lat))*SILx(:,pindx.lbb).*Si2C];

            Sixx(:,kk) = mfactor(FD, tmp);
            kk = kk + 1;
        end
        % bt bt
        if (par.opt_bt == on)
            [~,~,Hout] = buildPFD(par, 'bSi');
            PFD_bt_bt = Hout.PFD_bt_bt;
            k_bt_bt = d0(kappa_si./T.^2);
            par.k_bt_bt = k_bt_bt;
            par.PFD_bt_bt = PFD_bt_bt;
            tmp = bt*[k_bt*DSI; ...
                      -(k_bt+PFD_bt)*DSI] + ...
                  bt*bt*[k_bt_bt*DSI; ...
                         -(k_bt_bt+PFD_bt_bt)*DSI] + ...
                  2*bt*[k_bt*DSIx(:,pindx.lbt); ...
                        -(k_bt+PFD_bt)*DSIx(:,pindx.lbt)];
            
            Sixx(:,kk) = mfactor(FD, tmp);
            kk = kk + 1;
        end
        % bt aa
        if (par.opt_bt == on & par.opt_aa == on)
            tmp = bt*[k_bt*DSIx(:,pindx.aa); ...
                      -(k_bt + PFD_bt)*DSIx(:,pindx.aa)] + ...
                  [-d0(G*SILx(:,pindx.lbt))*dSi2Cdaa; ... 
                   d0(G*SILx(:,pindx.lbt))*dSi2Cdaa] + ...
                  [-d0(Gx(:,pindx.lbt))*SILx(:,pindx.aa).*Si2C; ...
                   d0(Gx(:,pindx.lbt))*SILx(:,pindx.aa).*Si2C];
            
            Sixx(:,kk) = mfactor(FD, tmp);
            kk = kk + 1;
        end
        % bt bb
        if (par.opt_bt == on & par.opt_bb == on)
            tmp = bt*[k_bt*DSIx(:,pindx.lbb); ...
                      -(k_bt+PFD_bt)*DSIx(:,pindx.lbb)] + ...
                  bb*[-d0(G*SILx(:,pindx.lbt))*dSi2Cdbb; ... 
                      d0(G*SILx(:,pindx.lbt))*dSi2Cdbb] + ...
                  [-d0(Gx(:,pindx.lbt))*SILx(:,pindx.lbb).*Si2C; ...
                   d0(Gx(:,pindx.lbt))*SILx(:,pindx.lbb).*Si2C];

            Sixx(:,kk) = mfactor(FD, tmp);
            kk = kk + 1;
        end
        % aa aa
        if (par.opt_aa == on & par.opt_aa == on)
            tmp = [-d0(Gxx(:,kk))*SIL; ...
                   d0(Gxx(:,kk))*SIL] + ...
                  2*[-d0(G*SILx(:,pindx.aa))*dSi2Cdaa; ...
                     d0(G*SILx(:,pindx.aa))*dSi2Cdaa];
            
            Sixx(:,kk) = mfactor(FD, tmp);
            kk = kk + 1;
        end
        % aa bb
        if (par.opt_aa == on & par.opt_bb == on)
            tmp = [-d0(Gxx(:,kk))*SIL; ...
                   d0(Gxx(:,kk))*SIL] + ...
                  bb*[-d0(G*SILx(:,pindx.aa))*dSi2Cdbb; ... 
                      d0(G*SILx(:,pindx.aa))*dSi2Cdbb] + ...
                  [-d0(G*SILx(:,pindx.lbb))*dSi2Cdaa; ... 
                   d0(G*SILx(:,pindx.lbb))*dSi2Cdaa];
            
            Sixx(:,kk) = mfactor(FD, tmp);
            kk = kk + 1;
        end
        % bb bb
        if (par.opt_bb == on & par.opt_bb == on)
            tmp = [-d0(Gxx(:,kk))*SIL; ...
                   d0(Gxx(:,kk))*SIL] + ...
                  2*bb*[-d0(G*SILx(:,pindx.lbb))*dSi2Cdbb; ...
                        d0(G*SILx(:,pindx.lbb))*dSi2Cdbb];
            
            Sixx(:,kk) = mfactor(FD, tmp);
            kk = kk + 1;
        end
    end
end

