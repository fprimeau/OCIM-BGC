function [par,P,Px,Pxx] = eqPcycle(x, par)
% ip is the mapping from x to parameter names (see switch below)
% output: P is model prediction of DIP,POP,and DOP
% output: F partial derivative of P model w.r.t. model parameters x
% output: Fxx hessian matrix of P model w.r.t.  model parameters x
% unpack some useful stuff
    on = true; off = false;
    TRdiv = par.TRdiv ;
    M3d   = par.M3d   ;
    grd   = par.grd   ;
    pindx = par.pindx ;
    iwet  = par.iwet  ;
    nwet  = par.nwet  ; % number of wet points;
    I     = par.I     ; % make an identity matrix;

    % unpack the parameters to be optimized
    if (par.opt_sigma == on)
        lsigma = x(pindx.lsigma);
        sigma  = exp(lsigma);
    else
        sigma  = par.sigma;
    end
    par.sigma = sigma; % pass parameter to C/Si/O models

    %
    if (par.opt_kP_T == on)
        kP_T = x(pindx.kP_T);
    else
        kP_T = par.kP_T;
    end
    par.kP_T = kP_T; % pass parameter to C/Si/O models

    %
    if (par.opt_kdP == on)
        lkdP = x(pindx.lkdP);
        kdP  = exp(lkdP);
    else
        kdP  = par.kdP;
    end
    par.kdP  = kdP; % pass parameter to C/Si/O models

    %
    if (par.opt_bP_T == on)
        bP_T = x(pindx.bP_T);
    else
        bP_T = par.bP_T;
    end
    par.bP_T = bP_T; % pass parameter to C/Si/O models

    %
    if (par.opt_bP == on)
        lbP = x(pindx.lbP);
        bP  = exp(lbP);
    else
        bP  = par.bP;
    end
    par.bP  = bP; % pass parameter to C/Si/O models

    %
    if (par.opt_alpha == on)
        lalpha = x(pindx.lalpha);
        alpha  = exp(lalpha);
    else
        alpha  = par.alpha;
    end
    par.alpha  = alpha; % pass parameter to C/Si/O models

    %
    if (par.opt_beta == on)
        lbeta = x(pindx.lbeta);
        beta  = exp(lbeta);
    else
        beta  = par.beta;
    end
    par.beta  = beta; % pass parameter to C/Si/O models

    % fixed parameters
    Tz      = par.Tz ;
    DIPbar  = M3d(iwet)*par.DIPbar;  % gobal average PO4 conc.[umol/kg] [need to multiply by par.permil to get mmol m^-3];
    kappa_g = par.kappa_g; % PO4 geological restore const.[s^-1];
    kappa_p = par.kappa_p; % POP solubilization rate constant
    npp     = par.npp;     % net primary production
    npp1    = par.npp1 ;
    npp2    = par.npp2 ;
    kP      = kP_T*Tz + kdP ;
    % build part of the biological DIP uptake operator
    Lambda     = par.Lambda;
    LAM        = 0*M3d;
    LAM(:,:,1) = (npp1.^beta).*Lambda(:,:,1);
    LAM(:,:,2) = (npp2.^beta).*Lambda(:,:,2);
    L          = d0(LAM(iwet));  % PO4 assimilation rate [s^-1];
    par.L      = L;
    PO4 = par.po4obs(iwet) ; % [umol/kg];
    % particle flux
    PFD = buildPFD(par,'POP');
    junk = M3d ;
    junk(:,:,2:end) = 0 ;
    isrf = find(junk(iwet)) ;

    % build Jacobian equations.
    % column 1 dF/dDIP
    Fp{1,1} = TRdiv + alpha*L + kappa_g*I ;
    Fp{2,1} = -(1-sigma)*alpha*L;
    Fp{3,1} = -sigma*alpha*L;

    % column 2 dF/dPOP
    Fp{1,2} = 0*I;
    Fp{2,2} = PFD + kappa_p*I;
    Fp{3,2} = -kappa_p*I;

    % column 3 dF/dDOP
    Fp{1,3} = -d0(kP) ;
    Fp{2,3} = 0*I;
    Fp{3,3} = TRdiv + d0(kP) ;

    % right hand side of phosphate equations
    RHS = [kappa_g*DIPbar; ...
           sparse(nwet,1);...
           sparse(nwet,1)];

    fprintf('Solving P model ...\n') ;
	%keyboard
    % dP/dt + Fp*P = RHS
    % factoring Jacobian matrix
    FFp = mfactor(cell2mat(Fp));
    % solve for P-cycle model state
    P = mfactor(FFp,RHS);

    if (par.optim == off) | ~any([par.opt_sigma par.opt_kP_T par.opt_kdP par.opt_bP_T par.opt_bP par.opt_alpha par.opt_beta])
        Px = [];
    elseif (par.optim & nargout > 2)
        %
        % Compute the gradient of the solution wrt the parameters
        %
        Z   = sparse(nwet,1);
        DIP = P(1:nwet);
        POP = P(nwet+1:2*nwet);
        DOP = P(2*nwet+1:end);

        % sigma
        if (par.opt_sigma == on)
            tmp = sigma*[Z; alpha*L*DIP; -alpha*L*DIP];
            Px(:,pindx.lsigma) = mfactor(FFp,-tmp);
        end

        % kP_T
        if (par.opt_kP_T == on)
            kP_kP_T = Tz ;
            tmp = [-d0(kP_kP_T)*DOP; ...
                   Z; ...
                   d0(kP_kP_T)*DOP];
            Px(:,pindx.kP_T) = mfactor(FFp,-tmp);
        end

        % kdP
        if (par.opt_kdP == on)
            kP_kdP = kdP ;
            tmp = kP_kdP*[-DOP; Z; DOP];
            Px(:,pindx.lkdP) = mfactor(FFp,-tmp);
        end

        % bP_T
        if (par.opt_bP_T == on)
            [~,Gout] = buildPFD(par,'POP');
            PFD_bm = Gout.PFD_bm;
            tmp =  [Z; PFD_bm*POP; Z];
            Px(:,pindx.bP_T) = mfactor(FFp, -tmp);
        end

        % bP
        if (par.opt_bP == on)
            [~,Gout] = buildPFD(par,'POP');
            PFD_bb = Gout.PFD_bb;
            tmp = bP*[Z; PFD_bb*POP;  Z];
            Px(:,pindx.lbP) = mfactor(FFp,-tmp);
        end

        % alpha
        if (par.opt_alpha == on)
            tmp = alpha*[L*DIP;...
                         -(1-sigma)*L*DIP;...
                         -sigma*L*DIP];
            Px(:,pindx.lalpha) = mfactor(FFp,-tmp);
        end

        %beta
        if (par.opt_beta == on)
            dLambdadbeta = 0*Lambda;
            dLambdadbeta(:,:,1) = log(npp1).*LAM(:,:,1);
            dLambdadbeta(:,:,2) = log(npp2).*LAM(:,:,2);
            iz = find(isinf(dLambdadbeta(:)));
            dLambdadbeta(iz) = 0;
            inan = find(isnan(dLambdadbeta(:)));
            dLambdadbeta(inan) = 0;
            dLdbeta = d0(dLambdadbeta(iwet));
            tmp = beta*[ alpha*dLdbeta*DIP;...
                         -(1-sigma)*alpha*dLdbeta*DIP;...
                         -sigma*alpha*dLdbeta*DIP];
            % will need L and dLdbeta for gradients of other
            % biogeochemical cycles
            par.L = L;
            par.dLdbeta = dLdbeta;
            Px(:,pindx.lbeta) = mfactor(FFp,-tmp);
        end
    end
    %% ---------------------------------------------------------
    if (par.optim == off) | ~any([par.opt_sigma par.opt_kP_T par.opt_kdP par.opt_bP_T par.opt_bP par.opt_alpha par.opt_beta])
        Pxx = [];
    elseif (par.optim & nargout > 3)
        % Compute the hessian of the solution wrt the parameters
        DIPx = Px(0*nwet+1:1*nwet,:);
        POPx = Px(1*nwet+1:2*nwet,:);
        DOPx = Px(2*nwet+1:end,:);
        % sigma sigma
        kk = 1;
        if (par.opt_sigma == on)
            tmp = [Z; ... % d2Jdsigma2 * DIP
                   sigma*alpha*L*DIP; ...
                   -sigma*alpha*L*DIP] + ...
                  2*[Z; ... % 2 * dJdigmadPdsigma
                     sigma*alpha*L*DIPx(:,pindx.lsigma); ...
                     -sigma*alpha*L*DIPx(:,pindx.lsigma)];

            Pxx(:,kk) = mfactor(FFp,-tmp);
            kk = kk + 1;
        end

        % sigma kP_T
        if (par.opt_sigma == on & par.opt_kP_T == on)
            tmp = [Z; Z; Z] + ... % d2Jdsigmadkappa * DIP
                  [Z; ... % dJdsigmadPdkappa
                   sigma*alpha*L*DIPx(:,pindx.kP_T);...
                   -sigma*alpha*L*DIPx(:,pindx.kP_T)] + ...
                  [-d0(kP_kP_T)*DOPx(:,pindx.lsigma); ... % J_kP_T*P_sigma
                   Z; ...
                   d0(kP_kP_T)*DOPx(:,pindx.lsigma)];

            Pxx(:,kk) = mfactor(FFp,-tmp);
            kk = kk + 1;
        end

        % sigma kdP
        if (par.opt_sigma == on & par.opt_kdP == on)
            tmp = [Z; Z; Z] + ... % d2Jdsigmadkappa * DIP
                  [Z; ... % dJdsigmadPdkappa
                   sigma*alpha*L*DIPx(:,pindx.lkdP);...
                   -sigma*alpha*L*DIPx(:,pindx.lkdP)] + ...
                  [-kdP*I*DOPx(:,pindx.lsigma); ... % dJdkappadPdsigma
                   Z; ...
                   kdP*I*DOPx(:,pindx.lsigma)];

            Pxx(:,kk) = mfactor(FFp,-tmp);
            kk = kk + 1;
        end

        % sigma bP_T
        if (par.opt_sigma == on & par.opt_bP_T == on)
            tmp = [Z; Z; Z] + ... % d2Jdsigmadslope * P
                  [Z; ... % dJdsigmadPdslope
                   sigma*alpha*L*DIPx(:,pindx.bP_T); ...
                   -sigma*alpha*L*DIPx(:,pindx.bP_T)] + ...
                  [Z;...  % dJdbP_TdPdsigma
                   PFD_bm*POPx(:,pindx.lsigma); ...
                   Z];

            Pxx(:,kk) = mfactor(FFp,-tmp);
            kk = kk + 1;
        end

        % sigma bP
        if (par.opt_sigma == on & par.opt_bP == on)
            tmp = [Z; Z; Z] + ... % d2Jdsigmadinterp * P
                  [Z; ... % dJdsigmadPdslope
                   sigma*alpha*L*DIPx(:,pindx.lbP); ...
                   -sigma*alpha*L*DIPx(:,pindx.lbP)] + ...
                  [Z; ... % dJdbPdPdsigma
                   bP*PFD_bb*POPx(:,pindx.lsigma); ...
                   Z];

            Pxx(:,kk) = mfactor(FFp,-tmp);
            kk = kk + 1;
        end

        % sigma alpha
        if (par.opt_sigma == on & par.opt_alpha == on)
            tmp = [Z; ... % d2Jdsigmadalpha
                   sigma*alpha*L*DIP;...
                   -sigma*alpha*L*DIP] + ...
                  [Z; ... % dJdsigmadPdalpha
                   sigma*alpha*L*DIPx(:,pindx.lalpha); ...
                   -sigma*alpha*L*DIPx(:,pindx.lalpha)] + ...
                  [alpha*L*DIPx(:,pindx.lsigma);... % dJdalphadsigma
                   -(1-sigma)*alpha*L*DIPx(:,pindx.lsigma)
                   -sigma*alpha*L*DIPx(:,pindx.lsigma)];

            Pxx(:,kk) = mfactor(FFp,-tmp);
            kk = kk + 1;
        end

        % sigma beta
        if (par.opt_sigma == on & par.opt_beta == on)
            tmp = [Z; ... % d2Jdsigmadbeta
                   sigma*alpha*beta*dLdbeta*DIP;...
                   -sigma*alpha*beta*dLdbeta*DIP] + ...
                  [Z; ... %dJdsigmadPdbeta
                   sigma*alpha*L*DIPx(:,pindx.lbeta); ...
                   -sigma*alpha*L*DIPx(:,pindx.lbeta)] + ...
                  [beta*alpha*dLdbeta*DIPx(:,pindx.lsigma);...  % dJdbetadPdsigma
                   -(1-sigma)*alpha*beta*dLdbeta*DIPx(:,pindx.lsigma);...
                   -sigma*alpha*beta*dLdbeta*DIPx(:,pindx.lsigma)];;

            Pxx(:,kk) = mfactor(FFp,-tmp);
            kk = kk + 1;
        end

        % kP_T kP_T
        if (par.opt_kP_T == on)
            tmp = 2*[-d0(kP_kP_T)*DOPx(:,pindx.kP_T); ... % dJdkappadPdkappa
                     Z;...
                     d0(kP_kP_T)*DOPx(:,pindx.kP_T)];

            Pxx(:,kk) = mfactor(FFp,-tmp);
            kk = kk + 1;
        end

        % kP_T kdP
        if (par.opt_kP_T == on & par.opt_kdP == on)
            tmp = [Z; Z; Z] + ... % d2JdkP_Tdkappa
                  [-d0(kP_kP_T)*DOPx(:,pindx.lkdP); ... % dJdkP_TdPdkappa
                   Z;...
                   d0(kP_kP_T)*DOPx(:,pindx.lkdP)] + ...
                  [-kP_kdP*DOPx(:,pindx.kP_T); ...  % dJdkappadPdkP_T
                   Z; ...
                   kP_kdP*DOPx(:,pindx.kP_T)];

            Pxx(:,kk) = mfactor(FFp,-tmp);
            kk = kk + 1;
        end

        % kP_T bP_T
        if (par.opt_kP_T == on & par.opt_bP_T == on)
            tmp = [Z; Z; Z] + ... % d2JdkP_Tdslope
                  [-d0(kP_kP_T)*DOPx(:,pindx.bP_T); ... % dJdkP_TdPdslope
                   Z;...
                   d0(kP_kP_T)*DOPx(:,pindx.bP_T)] + ...
                  [Z; ...  % dJdslopedPdkP_T
                   PFD_bm*POPx(:,pindx.kP_T); ...
                   Z];

            Pxx(:,kk) = mfactor(FFp,-tmp);
            kk = kk + 1;
        end

        % kP_T bP
        if (par.opt_kP_T == on & par.opt_bP == on)
            tmp = [Z; Z; Z] + ... % d2JdkP_Tdinterp
                  [-d0(kP_kP_T)*DOPx(:,pindx.lbP); ... % dJdkP_TdPdinterp
                   Z;...
                   d0(kP_kP_T)*DOPx(:,pindx.lbP)] + ...
                  [Z; ...  % dJdinterpdPdkP_T
                   bP*PFD_bb*POPx(:,pindx.kP_T); ...
                   Z];

            Pxx(:,kk) = mfactor(FFp,-tmp);
            kk = kk + 1;
        end

        % kP_T alpha
        if (par.opt_kP_T == on & par.opt_alpha == on)
            tmp = [Z; Z; Z] + ... % d2JdkP_Tdalpha
                  [-d0(kP_kP_T)*DOPx(:,pindx.lalpha); ... % dJdkP_TdPdalpha
                   Z;...
                   d0(kP_kP_T)*DOPx(:,pindx.lalpha)] + ...
                  alpha*[L*DIPx(:,pindx.kP_T);...  % dJdalphadPdkP_T
                         -(1-sigma)*L*DIPx(:,pindx.kP_T);...
                         -sigma*L*DIPx(:,pindx.kP_T)];

            Pxx(:,kk) = mfactor(FFp,-tmp);
            kk = kk + 1;
        end

        % kP_T beta
        if (par.opt_kP_T == on & par.opt_beta == on)
            tmp = [Z; Z; Z] + ... % d2JdkP_Tdbeta
                  [-d0(kP_kP_T)*DOPx(:,pindx.lbeta); ... % dJdkP_TdPdbeta
                   Z;...
                   d0(kP_kP_T)*DOPx(:,pindx.lbeta)] + ...
                  beta*[alpha*dLdbeta*DIPx(:,pindx.kP_T);... % dJdbetadPdkP_T
                        -(1-sigma)*alpha*dLdbeta*DIPx(:,pindx.kP_T);...
                        -sigma*alpha*dLdbeta*DIPx(:,pindx.kP_T)];

            Pxx(:,kk) = mfactor(FFp,-tmp);
            kk = kk + 1;
        end

        % kdP kdP
        if (par.opt_kdP == on)
            tmp = [-kP_kdP*DOP; ... % d2Jdkappa2 * P
                   Z;...
                   kP_kdP*DOP] + ...
                  2*[-kP_kdP*DOPx(:,pindx.lkdP); ... % dJdkappadPdkappa
                     Z;...
                     kP_kdP*DOPx(:,pindx.lkdP)];

            Pxx(:,kk) = mfactor(FFp,-tmp);
            kk = kk + 1;
        end

        % kdP bP_T
        if (par.opt_kdP == on & par.opt_bP_T == on)
            tmp = [Z; Z; Z] + ... % d2Jdkappadslope
                  [-kP_kdP*DOPx(:,pindx.bP_T); ... % dJdkappadPdslope
                   Z;...
                   kP_kdP*DOPx(:,pindx.bP_T)] + ...
                  [Z; ...  % dJdslopedPdkappa
                   PFD_bm*POPx(:,pindx.lkdP); ...
                   Z];

            Pxx(:,kk) = mfactor(FFp,-tmp);
            kk = kk + 1;
        end

        % kdP bP
        if (par.opt_kdP == on & par.opt_bP == on)
            tmp = [Z; Z; Z] + ... % d2Jdkappadinterp
                  [-kP_kdP*DOPx(:,pindx.lbP); ... % dJdkappadPdinterp
                   Z;...
                   kP_kdP*DOPx(:,pindx.lbP)] + ...
                  [Z; ...  % dJdinterpdPdkappa
                   bP*PFD_bb*POPx(:,pindx.lkdP); ...
                   Z];

            Pxx(:,kk) = mfactor(FFp,-tmp);
            kk = kk + 1;
        end

        % kdP alpha
        if (par.opt_kdP == on & par.opt_alpha == on)
            tmp = [Z; Z; Z] + ... % d2Jdkappadalpha
                  kP_kdP*[-DOPx(:,pindx.lalpha); ... % dJdkappadPdalpha
                          Z;...
                          DOPx(:,pindx.lalpha)] + ...
                  alpha*[L*DIPx(:,pindx.lkdP);...  % dJdalphadPdkappa
                         -(1-sigma)*L*DIPx(:,pindx.lkdP);...
                         -sigma*L*DIPx(:,pindx.lkdP)];

            Pxx(:,kk) = mfactor(FFp,-tmp);
            kk = kk + 1;
        end

        % kdP beta
        if (par.opt_kdP == on & par.opt_beta == on)
            tmp = [Z; Z; Z] + ... % d2Jdkappadbeta
                  kP_kdP*[-DOPx(:,pindx.lbeta); ... % dJdkappadPdbeta
                          Z;...
                          DOPx(:,pindx.lbeta)] + ...
                  beta*[alpha*dLdbeta*DIPx(:,pindx.lkdP);... % dJdbetadPdkappa
                        -(1-sigma)*alpha*dLdbeta*DIPx(:,pindx.lkdP);...
                        -sigma*alpha*dLdbeta*DIPx(:,pindx.lkdP)];

            Pxx(:,kk) = mfactor(FFp,-tmp);
            kk = kk + 1;
        end

        % bP_T bP_T
        if (par.opt_bP_T == on)
            [~,~,Hout] = buildPFD(par,'POP');
            PFD_bm_bm = Hout.PFD_bm_bm;
            tmp = [Z; ... % d2Jdslope2
                   PFD_bm_bm*POP; ...
                   Z] + ...
                  2*[Z; ... % dJdslopedPslope
                     PFD_bm*POPx(:,pindx.bP_T);...
                     Z];

            Pxx(:,kk) = mfactor(FFp,-tmp);
            kk = kk + 1;
        end

        % bP_T bP
        if (par.opt_bP_T == on & par.opt_bP == on)
            [~,~,Hout] = buildPFD(par,'POP');
            PFD_bm_bb = Hout.PFD_bm_bb;
            tmp = bP*[Z; ...  % d2Jdslopedinterp
                      PFD_bm_bb*POP; ...
                      Z] + ...
                  [Z; ... % dJdslopedPdinterp
                   PFD_bm*POPx(:,pindx.lbP);...
                   Z] + ...
                  bP*[Z; ...  % dJdinterpdPdslope
                      PFD_bb*POPx(:,pindx.bP_T); ...
                      Z];

            Pxx(:,kk) = mfactor(FFp,-tmp);
            kk = kk + 1;
        end

        % bP_T alpha
        if (par.opt_bP_T == on & par.opt_alpha == on)
            tmp = [Z; Z; Z] + ... % d2Jdslopedalpha
                  [Z; ... % dJdslopedPdalpha
                   PFD_bm*POPx(:,pindx.lalpha);...
                   Z] + ...
                  [alpha*L*DIPx(:,pindx.bP_T); ...  % dJdalphadPdslope
                   -(1-sigma)*alpha*L*DIPx(:,pindx.bP_T); ...
                   -sigma*alpha*L*DIPx(:,pindx.bP_T)];

            Pxx(:,kk) = mfactor(FFp,-tmp);
            kk = kk + 1;
        end

        % bP_T beta
        if (par.opt_bP_T == on & par.opt_beta == on)
            tmp = [Z; Z; Z] + ... % d2Jdslopedbeta
                  [Z; ... % dJdslopedPdbeta
                   PFD_bm*POPx(:,pindx.lbeta);...
                   Z] + ...
                  [alpha*beta*dLdbeta*DIPx(:,pindx.bP_T);...  % dJdbetadPdkappa
                   -(1-sigma)*alpha*beta*dLdbeta*DIPx(:,pindx.bP_T);...
                   -sigma*alpha*beta*dLdbeta*DIPx(:,pindx.bP_T)];

            Pxx(:,kk) = mfactor(FFp,-tmp);
            kk = kk + 1;
        end

        % bP bP
        if (par.opt_bP == on)
            [~,~,Hout] = buildPFD(par,'POP');
            PFD_bb_bb = Hout.PFD_bb_bb;
            tmp = [Z; ... % d2Jdinterp2
                   bP*bP*PFD_bb_bb*POP; ...
                   Z] + ...
                  [Z; ... % d2Jdinterp2
                   bP*PFD_bb*POP; ...
                   Z] + ...
                  2*[Z; ... % dJdinterpdPinterp
                     bP*PFD_bb*POPx(:,pindx.lbP);...
                     Z];

            Pxx(:,kk) = mfactor(FFp,-tmp);
            kk = kk + 1;
        end

        % bP alpha
        if (par.opt_bP == on & par.opt_alpha == on)
            tmp = [Z; Z; Z] + ... % d2Jdinterpdalpha
                  [Z; ... % dJdinterpdPdalpha
                   bP*PFD_bb*POPx(:,pindx.lalpha);...
                   Z] + ...
                  [alpha*L*DIPx(:,pindx.lbP); ... % dJdalphadPdinterp
                   -(1-sigma)*alpha*L*DIPx(:,pindx.lbP); ...
                   -sigma*alpha*L*DIPx(:,pindx.lbP)];

            Pxx(:,kk) = mfactor(FFp,-tmp);
            kk = kk + 1;
        end

        % bP beta
        if (par.opt_bP == on & par.opt_beta == on)
            tmp = [Z; ...
                   Z; ...
                   Z] + ... % d2Jdinterpdbeta
                  [Z; ... % dJdinterpdPdbeta
                   bP*PFD_bb*POPx(:,pindx.lbeta);...
                   Z] + ...
                  beta*[alpha*dLdbeta*DIPx(:,pindx.lbP); ... % dJdbetadPdinterp
                        -(1-sigma)*alpha*dLdbeta*DIPx(:,pindx.lbP); ...
                        -sigma*alpha*dLdbeta*DIPx(:,pindx.lbP)];

            Pxx(:,kk) = mfactor(FFp,-tmp);
            kk = kk + 1;
        end

        % alpha alpha
        if (par.opt_alpha == on)
            tmp = [alpha*L*DIP;... % d2Jdalpha2*P
                   -(1-sigma)*alpha*L*DIP;...
                   -alpha*sigma*L*DIP] + ...
                  2*[alpha*L*DIPx(:,pindx.lalpha); ...  % dJdalphadPdalpha
                     -(1-sigma)*alpha*L*DIPx(:,pindx.lalpha); ...
                     -sigma*alpha*L*DIPx(:,pindx.lalpha)];

            Pxx(:,kk) = mfactor(FFp,-tmp);
            kk = kk + 1;
        end

        % alpha beta
        if (par.opt_alpha == on & par.opt_beta == on)
            tmp = [alpha*beta*dLdbeta*DIP;... % d2Jdbeta2*P
                   -(1-sigma)*alpha*beta*dLdbeta*DIP;...
                   -sigma*alpha*beta*dLdbeta*DIP] + ...
                  alpha*[L*DIPx(:,pindx.lbeta);... % dJdalphadPdbeta
                         -(1-sigma)*L*DIPx(:,pindx.lbeta);...
                         -sigma*L*DIPx(:,pindx.lbeta)] + ...
                  beta*[alpha*dLdbeta*DIPx(:,pindx.lalpha);... %dJdbetadPdalpha
                        -(1-sigma)*alpha*dLdbeta*DIPx(:,pindx.lalpha);...
                        -sigma*alpha*dLdbeta*DIPx(:,pindx.lalpha)];

            Pxx(:,kk) = mfactor(FFp,-tmp);
            kk = kk + 1;
        end

        % beta beta
        if (par.opt_beta == on)
            d2Lambdadbetadbeta = 0*Lambda;
            d2Lambdadbetadbeta(:,:,1) = log(npp1).*log(npp1).*LAM(:,:,1);
            d2Lambdadbetadbeta(:,:,2) = log(npp2).*log(npp2).*LAM(:,:,2);
            iz = find(isinf(d2Lambdadbetadbeta(:)));
            d2Lambdadbetadbeta(iz) = 0;
            inan = find(isnan(d2Lambdadbetadbeta(:)));
            d2Lambdadbetadbeta(inan) = 0;
            d2Ldbetadbeta = d0(d2Lambdadbetadbeta(iwet));
            par.d2Ldbetadbeta = d2Ldbetadbeta;
            tmp = [alpha*beta*dLdbeta*DIP;... % d2Jdbeta2 * P
                   -(1-sigma)*alpha*beta*dLdbeta*DIP;...
                   -sigma*alpha*beta*dLdbeta*DIP] + ...
                  [alpha*beta*beta*d2Ldbetadbeta*DIP;... % d2Jdbeta2 * P
                   -(1-sigma)*alpha*beta*beta*d2Ldbetadbeta*DIP; ...
                   -sigma*alpha*beta*beta*d2Ldbetadbeta*DIP] + ...
                  2*[ alpha*beta*dLdbeta*DIPx(:,pindx.lbeta);...
                      -(1-sigma)*alpha*beta*dLdbeta*DIPx(:,pindx.lbeta);...
                      -sigma*alpha*beta*dLdbeta*DIPx(:,pindx.lbeta)];

            Pxx(:,kk) = mfactor(FFp,-tmp);
            kk = kk + 1;
        end
    end
end
