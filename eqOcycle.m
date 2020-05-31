function [O2, Ox] = eqOcycle(par, parm, x)
    on = true; off = false;
    global GO
    % unpack the parameters to be optimized
    
    % kappa_dc
    if (par.biogeochem.opt_kappa_dc == on)
        lkappa_dc = x(par.pindx.lkappa_dc);
        parm.kappa_dc  = exp(lkappa_dc);
    else
        parm.kappa_dc  = par.biogeochem.kappa_dc;
    end

    % slopeo
    if (par.biogeochem.opt_slopeo == on)
        parm.slopeo = x(par.pindx.slopeo);
    else
        parm.slopeo  = par.biogeochem.slopeo;
    end
    
    % interpo
    if (par.biogeochem.opt_interpo == on)
        linterpo = x(par.pindx.linterpo);
        parm.interpo  = exp(linterpo);
    else
        parm.interpo  = par.biogeochem.interpo;
    end
    
    X0  = GO;
    options.iprint = 1;
    options.atol = 1e-8; options.rtol = 1e-8 ;
    [O2,ierr] = nsnew(X0,@(X) O_eqn(X,par,parm,x),options) ;

    if (ierr ~=0)
        fprintf('o2model did not converge.\n') ;
        keyboard
    else
        % reset the global variable for the next call eqCcycle
        GO = real(O2) + 1e-3*randn(parm.nwet,1);
        X0 = GO;
        [O2,ierr] = nsnew(X0,@(X) O_eqn(X,par,parm,x),options);
        if nargout>1     
            %
            % Compute the gradient of the solution wrt the parameters
            [F, FD, Ox] = O_eqn(O2, par, parm, x);
            %
        end
    end

function [F, FD, Ox] = O_eqn(O2, par, parm, x)
    on = true; off = false;
    % fixed parameters
    grd  = parm.grd;
    M3d  = parm.M3d;
    iwet = parm.iwet;
    nwet = parm.nwet;
    sst  = parm.sst;
    sal  = parm.ss;
    TRdiv = parm.TRdiv;
    I = speye(nwet);
    
    % tunable parameters;
    kappa_dc = parm.kappa_dc;
    slopeo   = parm.slopeo;
    interpo  = parm.interpo;

    % variables from C model
    DIC = parm.DIC;
    DOC = parm.DOC;
    DOCx = parm.DOCx;
    
    smsk = M3d;
    smsk(:,:,2:end) = 0;
    isrf = find(smsk(iwet));
    dVs = parm.dVt(iwet(isrf));
    surface_mean = @(x) sum(x(isrf).*dVs)/sum(dVs);

    % compute the mean of the regressor variable
    Z = sst(iwet);
    mu = surface_mean(Z);
    Delta = sqrt(surface_mean((Z-mu).^2));

    % standardize the regressor variables
    ZR = (Z-mu)/Delta; parm.ZR = ZR;
    %
    O2C = slopeo*ZR + interpo;

    % O2 saturation concentration
    [KO2,o2sat] = Fsea2air_o2(parm);
    
    % organic C production rate
    [G,Gx] = uptake(parm, par);
    
    % rate of o2 production
    PO2 = O2C.*G;

    % parobolic function for o2 consumption
    R = 0.5 + 0.5*tanh(O2-10);
    dRdO = 0.5 - 0.5*tanh(O2-10).^2;
    
    % rate of o2 utilization
    % DOC(DOC<0) = 0;
    LO2 = O2C.*(kappa_dc*DOC.*R);
    dLdO = d0(O2C.*(kappa_dc*DOC.*dRdO));
    
    % O2 function
    F = TRdiv * O2 - PO2 + LO2 - KO2 * (o2sat - O2);
    FD = mfactor(TRdiv + dLdO + KO2);
    
    if (nargout > 2);
        % interpp
        if (par.biogeochem.opt_interpp == on)
            tmp = -O2C.*Gx(:,par.pindx.linterpp) + ...
                  kappa_dc*d0(O2C.*DOCx(:,par.pindx.linterpp))*R;
            Ox(:,par.pindx.linterpp) = mfactor(FD, -tmp);
        end
        
        % slopep
        if (par.biogeochem.opt_slopep == on)
            tmp = -O2C.*Gx(:,par.pindx.slopep) + ...
                  kappa_dc*O2C.*R.*DOCx(:,par.pindx.slopep);
            Ox(:,par.pindx.slopp) = mfactor(FD, -tmp);
        end
        
        % interpc
        if (par.biogeochem.opt_interpc == on)
            tmp =  kappa_dc*O2C.*R.*DOCx(:,par.pindx.linterpc);
            Ox(:,par.pindx.linterpc) = mfactor(FD, -tmp);
        end
        
        % slopec
        if (par.biogeochem.opt_slopec == on)
            tmp = kappa_dc*O2C.*R.*DOCx(:,par.pindx.slopec);
            Ox(:,par.pindx.slopec) = mfactor(FD, -tmp);
        end
        
        % sigma
        if (par.biogeochem.opt_sigma == on)
            tmp = -O2C.*Gx(:,par.pindx.lsigma) + ...
                  kappa_dc*O2C.*R.*DOCx(:,par.pindx.lsigma);
            Ox(:,par.pindx.lsigma) = mfactor(FD, -tmp);
        end
        
        % kappa_dp
        if (par.biogeochem.opt_kappa_dp == on)
            tmp = -O2C.*Gx(:,par.pindx.lkappa_dp) + ...
                  kappa_dc*O2C.*R.*DOCx(:,par.pindx.lkappa_dp);
            Ox(:,par.pindx.lkappa_dp) = mfactor(FD, -tmp);
        end
        
        % alpha
        if (par.biogeochem.opt_alpha == on)
            tmp = -O2C.*Gx(:,par.pindx.lalpha) + ...
                  kappa_dc*O2C.*R.*DOCx(:,par.pindx.lalpha);
            Ox(:,par.pindx.lalpha) = mfactor(FD, -tmp);
        end
        
        % beta
        if (par.biogeochem.opt_beta == on)
            tmp = -O2C.*Gx(:,par.pindx.lbeta) + ...
                  kappa_dc*O2C.*R.*DOCx(:,par.pindx.lbeta);
            Ox(:,par.pindx.lbeta) = mfactor(FD, -tmp);
        end
        
        % d
        if (par.biogeochem.opt_d == on)
            tmp = kappa_dc*O2C.*R.*DOCx(:,par.pindx.ld);
            Ox(:,par.pindx.ld) = mfactor(FD, -tmp);
        end
        
        % kappa_dc
        if (par.biogeochem.opt_kappa_dc == on)
            tmp = kappa_dc*DOC.*O2C.*R + ...
                  kappa_dc*O2C.*R.*DOCx(:,par.pindx.lkappa_dc);
            Ox(:,par.pindx.lkappa_dc) = mfactor(FD, -tmp);
        end
        
        % RR
        if (par.biogeochem.opt_RR == on)
            tmp = kappa_dc*O2C.*R.*DOCx(:,par.pindx.lRR);
            Ox(:,par.pindx.lRR) = mfactor(FD, -tmp);
        end

        if (par.biogeochem.opt_slopeo == on);
            ZR = parm.ZR;
            tmp = -ZR.*G + kappa_dc*ZR.*DOC.*R;
            Ox(:,par.pindx.slopeo) = mfactor(FD, -tmp);
        end

        if (par.biogeochem.opt_interpo == on);
            tmp = -interpo*G + interpo*kappa_dc*DOC.*R;
            Ox(:,par.pindx.linterpo) = mfactor(FD, -tmp);
        end
    end