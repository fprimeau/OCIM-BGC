function [par,P,Px,Pxx] = eqPcycle(x, par)
% ip is the mapping from x to parameter names (see switch below)
% output: P is model prediction of DIP,POP,and DOP
% output: F partial derivative of P model w.r.t. model parameters x
% output: Fxx hessian matrix of P model w.r.t.  model parameters x
% unpack some useful stuff
    on = true; off = false;
    TRdiv = par.TRdiv ; % addvection diffusion transporter;
    M3d   = par.M3d   ; % wet-dry mask;
    grd   = par.grd   ; % grid info.;
    pindx = par.pindx ; % indices of parameters; 
    iwet  = par.iwet  ; % indices of wet points;
    nwet  = par.nwet  ; % number of wet points;
    I     = par.I     ; % make an identity matrix;

    % unpack the parameters to be optimized
    if (par.opt_sigP == on)
        lsigP = x(pindx.lsigP);
        sigP  = exp(lsigP);
    else
        sigP  = par.sigP;
    end
    par.sigP = sigP; % pass parameter to C/Si/O models

    %
    if (par.opt_Q10P == on)
        lQ10P = x( pindx.lQ10P ) ;
         Q10P = exp( lQ10P ) ;
    else
        Q10P = par.Q10P;
    end
    par.Q10P = Q10P; % pass parameter to C/Si/O models

    %
    if (par.opt_kdP == on)
        lkdP = x(pindx.lkdP);
        kdP  = exp(lkdP);
    else
        kdP  = par.kdP;
    end
    par.kdP  = kdP; % pass parameter to C/Si/O models

    % bP_T
    if (par.opt_bP_T == on)
        par.bP_T = x(par.pindx.bP_T);
    end

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
    vT      = par.vT ;
    DIPbar  = M3d(iwet)*par.DIPbar ;  % gobal arerage PO4 conc.[mmol m^-3];
    gamma   = par.gamma   ; % fraction goes to lDOM ;
    kappa_g = par.kappa_g ; % PO4 geological restore const.[s^-1] ;
    kappa_p = par.kappa_p ; % POP solubilization rate constant ;
    kappa_l = par.kappa_l ; % labile DOM remineralization rate [s^-1];

    npp     = par.npp ;     % net primary production
    tf      = (vT - 30)/10 ;
    kP      =  kdP * Q10P.^tf ;
    % build part of the biological DIP uptake operator
    Lambda = par.Lambda;
    LAM    = 0*M3d;

    % change NaN value to zero in npp and Lambda
    % ---> it happens due to DIP_obs has nan values
    % in Weilei's code, smoothit functon remove the nan values of PO4.

    Lambda(isnan(Lambda(:))) = 0;
    npp(isnan(npp(:))) = 0;

    for ji = 1 : par.nl
        LAM(:,:,ji) = (npp(:,:,ji).^beta).*Lambda(:,:,ji) ;
    end 

    L      = d0(LAM(iwet));  % PO4 assimilation rate [s^-1];
    par.L  = L;
   
    % particle flux
    PFD  = buildPFD(par,'POP');
    junk = M3d ;
    junk(:,:,2:end) = 0 ;
    isrf = find(junk(iwet)) ;
    
    % build Jacobian equations.
    % column 1 dF/dDIP
    Fp{1,1} = TRdiv + alpha*L + kappa_g*I ;
    Fp{2,1} = -(1-sigP-gamma)*alpha*L ;
    Fp{3,1} = -sigP*alpha*L  ;
    Fp{4,1} = -gamma*alpha*L  ;

    % column 2 dF/dPOP
    Fp{1,2} = -kappa_p*I      ;
    Fp{2,2} = PFD + kappa_p*I ;
    Fp{3,2} = 0*I             ;
    Fp{4,2} = 0*I             ;

    % column 3 dF/dDOP
    Fp{1,3} = -d0(kP)         ;
    Fp{2,3} = 0*I ;
    Fp{3,3} = TRdiv + d0(kP)  ;
    Fp{4,3} = 0*I ;

    % column 4 dF/dlDOP
    Fp{1,4} = -kappa_l*I ;
    Fp{2,4} = 0*I ;
    Fp{3,4} = 0*I;
    Fp{4,4} = TRdiv + kappa_l*I ;

    % right hand side of phosphate equations
    RHS = [kappa_g*DIPbar  ; ... 
           sparse(nwet,1)  ; ...
           sparse(nwet,1)  ; ...
           sparse(nwet,1)] ;

    fprintf('Solving P model ...\n') ;
    % dP/dt + Fp*P = RHS
    Fp = cell2mat(Fp) ;
    inan = find(isnan(Fp(:))) ;
    iinf = find(isinf(Fp(:))) ;
    if length(iinf) + length(inan) > 0
        fprintf('the LHS has nans or infs, please check')
    end
 
    fprintf('factoring Jacobian matrix ...\n')
    tic
    FFp = mfactor(Fp) ;
    toc

    fprintf('solve for P-cycle model state ...\n')
    tic
    P = mfactor(FFp, RHS) ;
    toc

    if (par.optim == off)
        Px = [];
    elseif (par.optim & nargout > 2)
        %
        fprintf('Compute the gradient of the solution wrt the parameters...\n')
        tic
        %
        Z   = sparse(nwet,1);
        DIP = P(0*nwet+1:1*nwet);
        POP = P(1*nwet+1:2*nwet);
        DOP = P(2*nwet+1:3*nwet);
        DOPl= P(3*nwet+1:4*nwet) ;
        % sigP
        if (par.opt_sigP == on)
            tmp = sigP*[Z; alpha*L*DIP; -alpha*L*DIP; Z];
            RHS(:,pindx.lsigP) = -tmp;
        end
        
        % Q10P
        if (par.opt_Q10P == on)
            kP_Q10P = kdP * Q10P .* Q10P.^(tf - 1) .* tf ; 
            tmp = [-d0(kP_Q10P)*DOP; ...
                   Z; ...
                   d0(kP_Q10P)*DOP; ...
                   Z];
            RHS(:,pindx.lQ10P) = -tmp;
        end
        
        % kdP
        if (par.opt_kdP == on)
            kP_kdP = kdP * d0(Q10P.^tf) ;
            tmp = [-kP_kdP*DOP; Z; kP_kdP*DOP; Z];
            RHS(:,pindx.lkdP) = -tmp;
        end

        % bP_T
        if (par.opt_bP_T == on)
            [~,Gout] = buildPFD_48layer(par,'POP');
            PFD_bm = Gout.PFD_bm;
            tmp =  [Z; PFD_bm*POP; Z; Z];
            RHS(:,pindx.bP_T) = -tmp;
        end

        % bP
        if (par.opt_bP == on)
            [~,Gout] = buildPFD_48layer(par,'POP');
            PFD_bb = Gout.PFD_bb;
            tmp = bP*[Z; PFD_bb*POP;  Z; Z];
            RHS(:,pindx.lbP) = -tmp;
        end
        
        % alpha
        if (par.opt_alpha == on)
            tmp = alpha*[L*DIP; ...
                         -(1-sigP-gamma)*L*DIP; ...
                         -sigP*L*DIP; ...
                         -gamma*L*DIP];
            RHS(:,pindx.lalpha) = -tmp;
        end
        
        %beta
        if (par.opt_beta == on)
            dLambdadbeta = 0*Lambda;
            for ji = 1 : par.nl
                dLambdadbeta(:,:,ji) = log(npp(:,:,ji)).*LAM(:,:,ji) ;
            end 
            
            iz = find(isinf(dLambdadbeta(:)));
            dLambdadbeta(iz) = 0;
            inan = find(isnan(dLambdadbeta(:)));
            dLambdadbeta(inan) = 0;
            dLdbeta = d0(dLambdadbeta(iwet));
            tmp = beta*[ alpha*dLdbeta*DIP ; ...
                         -(1-sigP-gamma)*alpha*dLdbeta*DIP ; ...
                         -sigP*alpha*dLdbeta*DIP ; ...
                         -gamma*alpha*dLdbeta*DIP ];
            % will need L and dLdbeta for gradients of other
            % biogeochemical cycles
            par.L = L;
            par.dLdbeta = dLdbeta;
            RHS(:,pindx.lbeta) = -tmp;
        end
        Px = mfactor(FFp, RHS);
        toc
    end
    %% ---------------------------------------------------------
    if (par.optim == off)
        Pxx = [];
    elseif (par.optim & nargout > 3) 
        fprintf(' Compute the hessian of the solution wrt the parameters...\n')
        tic
        DIPx  = Px(0*nwet+1:1*nwet, :) ;
        POPx  = Px(1*nwet+1:2*nwet, :) ;
        DOPx  = Px(2*nwet+1:3*nwet, :) ;
        DOPlx = Px(3*nwet+1:4*nwet, :) ;
        % sigP sigP
        kk = 1;
        if (par.opt_sigP == on)
            tmp = [Z; ... % d2JdsigP2 * DIP 
                   sigP*alpha*L*DIP; ...
                   -sigP*alpha*L*DIP; ...
                   Z] + ...
                  2*[Z; ... % 2 * dJdigmadPdsigP
                     sigP*alpha*L*DIPx(:,pindx.lsigP); ...
                     -sigP*alpha*L*DIPx(:,pindx.lsigP);
                     Z];

            RHS(:,kk) = -tmp;
            kk = kk + 1;
        end

        % sigP Q10P
        if (par.opt_sigP == on & par.opt_Q10P == on)
            tmp = [Z; ... % dJdsigPdPdkappa
                   sigP*alpha*L*DIPx(:,pindx.lQ10P); ...
                   -sigP*alpha*L*DIPx(:,pindx.lQ10P); ...
                   Z] + ...
                  [-d0(kP_Q10P)*DOPx(:,pindx.lsigP); ... % J_Q10P*P_sigP
                   Z; ...
                   d0(kP_Q10P)*DOPx(:,pindx.lsigP); ...
                   Z];

            RHS(:,kk) = -tmp;
            kk = kk + 1;
        end

        % sigP kdP
        if (par.opt_sigP == on & par.opt_kdP == on)
            tmp = [Z; ... % dJdsigPdPdkappa
                   sigP*alpha*L*DIPx(:,pindx.lkdP); ...
                   -sigP*alpha*L*DIPx(:,pindx.lkdP); ...
                   Z] + ...
                  [-kP_kdP*DOPx(:,pindx.lsigP); ... % dJdkappadPdsigP
                   Z; ...
                   kP_kdP*DOPx(:,pindx.lsigP); ...
                   Z];

            RHS(:,kk) = -tmp;
            kk = kk + 1;
        end

        % sigP bP_T
        if (par.opt_sigP == on & par.opt_bP_T == on)
            tmp = [Z; ... % dJdsigPdPdslope
                   sigP*alpha*L*DIPx(:,pindx.bP_T); ...
                   -sigP*alpha*L*DIPx(:,pindx.bP_T); ...
                   Z] + ...
                  [Z;...  % dJdbP_TdPdsigP
                   PFD_bm*POPx(:,pindx.lsigP); ...
                   Z;...
                   Z];
            
            RHS(:,kk) = -tmp;
            kk = kk + 1;
        end

        % sigP bP
        if (par.opt_sigP == on & par.opt_bP == on)
            tmp = [Z; ... % dJdsigPdPdslope
                   sigP*alpha*L*DIPx(:,pindx.lbP); ...
                   -sigP*alpha*L*DIPx(:,pindx.lbP); ...
                   Z] + ...
                  [Z; ... % dJdbPdPdsigP
                   bP*PFD_bb*POPx(:,pindx.lsigP); ...
                   Z; ...
                   Z]; 
            
            RHS(:,kk) = -tmp;
            kk = kk + 1;
        end

        % sigP alpha
        if (par.opt_sigP == on & par.opt_alpha == on)
            tmp = [Z; ... % d2JdsigPdalpha
                   sigP*alpha*L*DIP;...
                   -sigP*alpha*L*DIP; ...
                   Z] + ...
                  [Z; ... % dJdsigPdPdalpha
                   sigP*alpha*L*DIPx(:,pindx.lalpha); ...
                   -sigP*alpha*L*DIPx(:,pindx.lalpha); ...
                   Z] + ...
                  [alpha*L*DIPx(:,pindx.lsigP);... % dJdalphadsigP
                   -(1-sigP-gamma)*alpha*L*DIPx(:,pindx.lsigP)
                   -sigP*alpha*L*DIPx(:,pindx.lsigP) ;
                   -gamma*alpha*L*DIPx(:,pindx.lsigP)];
            
            RHS(:,kk) = -tmp;
            kk = kk + 1;
        end

        % sigP beta
        if (par.opt_sigP == on & par.opt_beta == on)
            tmp = [Z; ... % d2JdsigPdbeta
                   sigP*alpha*beta*dLdbeta*DIP;...
                   -sigP*alpha*beta*dLdbeta*DIP; ...
                  Z] + ...
                  [Z; ... %dJdsigPdPdbeta
                   sigP*alpha*L*DIPx(:,pindx.lbeta); ...
                   -sigP*alpha*L*DIPx(:,pindx.lbeta); ...
                   Z] + ...
                  [beta*alpha*dLdbeta*DIPx(:,pindx.lsigP);...  % dJdbetadPdsigP
                   -(1-sigP-gamma)*alpha*beta*dLdbeta*DIPx(:,pindx.lsigP);...
                   -sigP*alpha*beta*dLdbeta*DIPx(:,pindx.lsigP) ;
                   -gamma*alpha*beta*dLdbeta*DIPx(:,pindx.lsigP)];
            
            RHS(:,kk) = -tmp;
            kk = kk + 1;
        end
        
        % Q10P Q10P
        if (par.opt_Q10P == on)
            kP_Q10P_Q10P = kP_Q10P + tf.*Q10P.^2.*kdP.*Q10P.^(tf - 2).*(tf - 1) ;
            tmp = 2*[-d0(kP_Q10P)*DOPx(:,pindx.lQ10P); ... % dJdkappadPdkappa
                     Z;...
                     d0(kP_Q10P)*DOPx(:,pindx.lQ10P); ...
                     Z] + ...
                  [-d0(kP_Q10P_Q10P)*DOP; ...
                   Z; ...
                   d0(kP_Q10P_Q10P)*DOP; ...
                   Z];
            
            RHS(:,kk) = -tmp;
            kk = kk + 1;
        end

        % Q10P kdP
        if (par.opt_Q10P == on & par.opt_kdP == on)
            kP_kdP_Q10P = tf .* Q10P .* kdP .* Q10P.^(tf-1);
            tmp = [-d0(kP_Q10P)*DOPx(:,pindx.lkdP); ... % dJdQ10PdPdkappa
                   Z;...
                   d0(kP_Q10P)*DOPx(:,pindx.lkdP); ...
                   Z] + ...
                  [-kP_kdP*DOPx(:,pindx.lQ10P); ...  % dJdkappadPdQ10P
                   Z; ... 
                   kP_kdP*DOPx(:,pindx.lQ10P); ...
                   Z] + ...
                  [-d0(kP_kdP_Q10P)*DOP; ...
                   Z; ...
                   d0(kP_kdP_Q10P)*DOP ; ...
                   Z]; 
            
            RHS(:, kk) = -tmp;
            kk = kk + 1;
        end
        
        % Q10P bP_T
        if (par.opt_Q10P == on & par.opt_bP_T == on)
            tmp = [-d0(kP_Q10P)*DOPx(:,pindx.bP_T); ... % dJdQ10PdPdslope
                   Z;...
                   d0(kP_Q10P)*DOPx(:,pindx.bP_T) ; ...
                   Z] + ...
                  [Z; ...  % dJdslopedPdQ10P
                   PFD_bm*POPx(:,pindx.lQ10P); ...
                   Z; ...
                   Z]; 
            
            RHS(:,kk) = -tmp;
            kk = kk + 1;
        end

        % Q10P bP
        if (par.opt_Q10P == on & par.opt_bP == on)
            tmp = [-d0(kP_Q10P)*DOPx(:,pindx.lbP); ... % dJdQ10PdPdinterp
                   Z;...
                   d0(kP_Q10P)*DOPx(:,pindx.lbP); ...
                   Z] + ...
                  [Z; ...  % dJdinterpdPdQ10P
                   bP*PFD_bb*POPx(:,pindx.lQ10P); ...
                   Z; ...
                   Z]; 
            
            RHS(:,kk) = -tmp;
            kk = kk + 1;
        end

        % Q10P alpha
        if (par.opt_Q10P == on & par.opt_alpha == on)
            tmp = [-d0(kP_Q10P)*DOPx(:,pindx.lalpha); ... % dJdQ10PdPdalpha
                   Z;...
                   d0(kP_Q10P)*DOPx(:,pindx.lalpha); ...
                   Z] + ...
                  alpha*[L*DIPx(:,pindx.lQ10P);...  % dJdalphadPdQ10P
                         -(1-sigP-gamma)*L*DIPx(:,pindx.lQ10P);...
                         -sigP*L*DIPx(:,pindx.lQ10P); ...
                         -gamma*L*DIPx(:,pindx.lQ10P)];
            
            RHS(:,kk) = -tmp;
            kk = kk + 1;
        end

        % Q10P beta
        if (par.opt_Q10P == on & par.opt_beta == on)
            tmp = [-d0(kP_Q10P)*DOPx(:,pindx.lbeta); ... % dJdQ10PdPdbeta
                   Z;...
                   d0(kP_Q10P)*DOPx(:,pindx.lbeta); ...
                   Z] + ...
                  beta*[alpha*dLdbeta*DIPx(:,pindx.lQ10P);... % dJdbetadPdQ10P
                        -(1-sigP-gamma)*alpha*dLdbeta*DIPx(:,pindx.lQ10P);...
                        -sigP*alpha*dLdbeta*DIPx(:,pindx.lQ10P); ...
                        -gamma*alpha*dLdbeta*DIPx(:,pindx.lQ10P)];
            
            RHS(:,kk) = -tmp;
            kk = kk + 1;
        end
        
        % kdP kdP
        if (par.opt_kdP == on)
            kP_kdP_kdP = kP_kdP ;
            tmp = [-kP_kdP_kdP*DOP; ... % d2Jdkappa2 * P
                   Z;...
                   kP_kdP_kdP*DOP; ...
                   Z] + ...
                  2*[-kP_kdP*DOPx(:,pindx.lkdP); ... % dJdkappadPdkappa
                     Z;...
                     kP_kdP*DOPx(:,pindx.lkdP); ...
                     Z];
            
            RHS(:,kk) = -tmp;
            kk = kk + 1;
        end

        % kdP bP_T
        if (par.opt_kdP == on & par.opt_bP_T == on)
            tmp = [-kP_kdP*DOPx(:,pindx.bP_T); ... % dJdkappadPdslope
                   Z;...
                   kP_kdP*DOPx(:,pindx.bP_T); ...
                   Z] + ...
                  [Z; ...  % dJdslopedPdkappa
                   PFD_bm*POPx(:,pindx.lkdP); ...
                   Z; ...
                   Z]; 
            
            RHS(:,kk) = -tmp;
            kk = kk + 1;
        end

        % kdP bP
        if (par.opt_kdP == on & par.opt_bP == on)
            tmp = [-kP_kdP*DOPx(:,pindx.lbP); ... % dJdkappadPdinterp
                   Z;...
                   kP_kdP*DOPx(:,pindx.lbP); ...
                   Z] + ...
                  [Z; ...  % dJdinterpdPdkappa
                   bP*PFD_bb*POPx(:,pindx.lkdP); ...
                   Z; ...
                   Z]; 
            
            RHS(:,kk) = -tmp;
            kk = kk + 1;
        end

        % kdP alpha
        if (par.opt_kdP == on & par.opt_alpha == on)
            tmp = [-kP_kdP*DOPx(:,pindx.lalpha); ... % dJdkappadPdalpha
                   Z;...
                   kP_kdP*DOPx(:,pindx.lalpha); ...
                   Z] + ...
                  alpha*[L*DIPx(:,pindx.lkdP);...  % dJdalphadPdkappa
                         -(1-sigP-gamma)*L*DIPx(:,pindx.lkdP);...
                         -sigP*L*DIPx(:,pindx.lkdP); ...
                         -gamma*L*DIPx(:,pindx.lkdP)];
            
            RHS(:,kk) = -tmp;
            kk = kk + 1;
        end

        % kdP beta
        if (par.opt_kdP == on & par.opt_beta == on)
            tmp = [-kP_kdP*DOPx(:,pindx.lbeta); ... % dJdkappadPdbeta
                   Z;...
                   kP_kdP*DOPx(:,pindx.lbeta); ...
                   Z] + ...
                  beta*[alpha*dLdbeta*DIPx(:,pindx.lkdP);... % dJdbetadPdkappa
                        -(1-sigP-gamma)*alpha*dLdbeta*DIPx(:,pindx.lkdP);...
                        -sigP*alpha*dLdbeta*DIPx(:,pindx.lkdP); ...
                        -gamma*alpha*dLdbeta*DIPx(:,pindx.lkdP)];
            
            RHS(:,kk) = -tmp;
            kk = kk + 1;
        end

        % bP_T bP_T
        if (par.opt_bP_T == on)
            [~,~,Hout] = buildPFD_48layer(par,'POP');
            PFD_bm_bm = Hout.PFD_bm_bm;
            tmp = [Z; ...
                   PFD_bm_bm*POP; ...
                   Z; ...
                   Z] + ...% d2Jdslope2
                  [Z; ...
                   2*PFD_bm*POPx(:,pindx.bP_T); ...
                   Z; ...
                   Z];% dJdslopedPslope
            
            RHS(:,kk) = -tmp;
            kk = kk + 1;
        end
        
        % bP_T bP
        if (par.opt_bP_T == on & par.opt_bP == on)
            [~,~,Hout] = buildPFD_48layer(par,'POP');
            PFD_bm_bb = Hout.PFD_bm_bb;
            tmp = [Z; ...  % d2Jdslopedinterp
                   bP*PFD_bm_bb*POP; ...
                   Z; ...
                   Z] + ...
                  [Z; ... % dJdslopedPdinterp
                   PFD_bm*POPx(:,pindx.lbP);...
                   Z; ...
                   Z] + ...
                  [Z; ...  % dJdinterpdPdslope
                   bP*PFD_bb*POPx(:,pindx.bP_T); ...
                   Z; ...
                   Z]; 
            
            RHS(:,kk) = -tmp;
            kk = kk + 1;
        end
        
        % bP_T alpha
        if (par.opt_bP_T == on & par.opt_alpha == on)
            tmp = [Z; ... % dJdslopedPdalpha
                   PFD_bm*POPx(:,pindx.lalpha);...
                   Z; ...
                   Z] + ...
                  [alpha*L*DIPx(:,pindx.bP_T); ...  % dJdalphadPdslope
                   -(1-sigP-gamma)*alpha*L*DIPx(:,pindx.bP_T); ...
                   -sigP*alpha*L*DIPx(:,pindx.bP_T); ...
                   -gamma*alpha*L*DIPx(:,pindx.bP_T)];
            
            RHS(:,kk) = -tmp;
            kk = kk + 1;
        end
        
        % bP_T beta
        if (par.opt_bP_T == on & par.opt_beta == on)
            tmp = [Z; ... % dJdslopedPdbeta
                   PFD_bm*POPx(:,pindx.lbeta);...
                   Z; ...
                   Z] + ...
                  [alpha*beta*dLdbeta*DIPx(:,pindx.bP_T);...  % dJdbetadPdkappa
                   -(1-sigP-gamma)*alpha*beta*dLdbeta*DIPx(:,pindx.bP_T);...
                   -sigP*alpha*beta*dLdbeta*DIPx(:,pindx.bP_T); ...
                   -gamma*alpha*beta*dLdbeta*DIPx(:,pindx.bP_T)];
            
            RHS(:,kk) = -tmp;
            kk = kk + 1;
        end
        
        % bP bP
        if (par.opt_bP == on)
            [~,~,Hout] = buildPFD_48layer(par,'POP');
            PFD_bb_bb = Hout.PFD_bb_bb;
            tmp = [Z; ... % d2Jdinterp2
                   bP*bP*PFD_bb_bb*POP; ...
                   Z; ...
                   Z] + ...
                  [Z; ... % d2Jdinterp2
                   bP*PFD_bb*POP; ...
                   Z; ...
                   Z] + ...
                  [Z; ... % dJdinterpdPinterp
                   2*bP*PFD_bb*POPx(:,pindx.lbP);...
                   Z; ...
                   Z];
            
            RHS(:,kk) = -tmp;
            kk = kk + 1;
        end    

        % bP alpha
        if (par.opt_bP == on & par.opt_alpha == on)
            tmp = [Z; ... % dJdinterpdPdalpha
                   bP*PFD_bb*POPx(:,pindx.lalpha);...
                   Z; ...
                   Z] + ...
                  [alpha*L*DIPx(:,pindx.lbP); ... % dJdalphadPdinterp
                   -(1-sigP-gamma)*alpha*L*DIPx(:,pindx.lbP); ...
                   -sigP*alpha*L*DIPx(:,pindx.lbP); ...
                   -gamma*alpha*L*DIPx(:,pindx.lbP)];  
            
            RHS(:,kk) = -tmp;
            kk = kk + 1;
        end

        % bP beta
        if (par.opt_bP == on & par.opt_beta == on)
            tmp = [Z; ... % dJdinterpdPdbeta
                   bP*PFD_bb*POPx(:,pindx.lbeta);...
                   Z; ...
                   Z] + ...
                  beta*[alpha*dLdbeta*DIPx(:,pindx.lbP); ... % dJdbetadPdinterp
                        -(1-sigP-gamma)*alpha*dLdbeta*DIPx(:,pindx.lbP); ...
                        -sigP*alpha*dLdbeta*DIPx(:,pindx.lbP); ...
                        -gamma*alpha*dLdbeta*DIPx(:,pindx.lbP)];
            
            RHS(:,kk) = -tmp;
            kk = kk + 1;
        end
        
        % alpha alpha
        if (par.opt_alpha == on)
            tmp = [alpha*L*DIP;... % d2Jdalpha2*P
                   -(1-sigP-gamma)*alpha*L*DIP;...
                   -alpha*sigP*L*DIP; ...
                   -alpha*gamma*L*DIP] + ...
                  2*[alpha*L*DIPx(:,pindx.lalpha); ...  % dJdalphadPdalpha
                     -(1-sigP-gamma)*alpha*L*DIPx(:,pindx.lalpha); ...
                     -sigP*alpha*L*DIPx(:,pindx.lalpha); ...
                     -gamma*alpha*L*DIPx(:,pindx.lalpha)];

            RHS(:,kk) = -tmp;
            kk = kk + 1;
        end

        % alpha beta
        if (par.opt_alpha == on & par.opt_beta == on)
            tmp = [alpha*beta*dLdbeta*DIP;... % d2Jdbeta2*P
                   -(1-sigP-gamma)*alpha*beta*dLdbeta*DIP;...
                   -sigP*alpha*beta*dLdbeta*DIP; ...
                   -gamma*alpha*beta*dLdbeta*DIP] + ...
                  alpha*[L*DIPx(:,pindx.lbeta);... % dJdalphadPdbeta
                         -(1-sigP-gamma)*L*DIPx(:,pindx.lbeta);...
                         -sigP*L*DIPx(:,pindx.lbeta); ...
                         -gamma*L*DIPx(:,pindx.lbeta)] + ...
                  beta*[alpha*dLdbeta*DIPx(:,pindx.lalpha);... %dJdbetadPdalpha
                        -(1-sigP-gamma)*alpha*dLdbeta*DIPx(:,pindx.lalpha);...
                        -sigP*alpha*dLdbeta*DIPx(:,pindx.lalpha); ...
                        -gamma*alpha*dLdbeta*DIPx(:,pindx.lalpha)];

            RHS(:,kk) = -tmp;
            kk = kk + 1;
        end
        
        % beta beta
        if (par.opt_beta == on)
            d2Lambdadbetadbeta = 0*Lambda;
            for ji = 1 : par.nl
                d2Lambdadbetadbeta(:,:,ji) = log(npp(:,:,ji)).*log(npp(:,:,ji)).*LAM(:,:,ji);
            end 
            
            iz = find(isinf(d2Lambdadbetadbeta(:)));
            d2Lambdadbetadbeta(iz) = 0;
            inan = find(isnan(d2Lambdadbetadbeta(:)));
            d2Lambdadbetadbeta(inan) = 0;
            d2Ldbetadbeta = d0(d2Lambdadbetadbeta(iwet));
            par.d2Ldbetadbeta = d2Ldbetadbeta;
            tmp = [alpha*beta*dLdbeta*DIP;... % d2Jdbeta2 * P
                   -(1-sigP-gamma)*alpha*beta*dLdbeta*DIP;...
                   -sigP*alpha*beta*dLdbeta*DIP; ...
                   -gamma*alpha*beta*dLdbeta*DIP] + ...
                  [alpha*beta*beta*d2Ldbetadbeta*DIP;... % d2Jdbeta2 * P
                   -(1-sigP-gamma)*alpha*beta*beta*d2Ldbetadbeta*DIP; ...
                   -sigP*alpha*beta*beta*d2Ldbetadbeta*DIP; ...
                   -gamma*alpha*beta*beta*d2Ldbetadbeta*DIP] + ...
                  2*[ alpha*beta*dLdbeta*DIPx(:,pindx.lbeta);...
                      -(1-sigP-gamma)*alpha*beta*dLdbeta*DIPx(:,pindx.lbeta);...
                      -sigP*alpha*beta*dLdbeta*DIPx(:,pindx.lbeta); ...
                      -gamma*alpha*beta*dLdbeta*DIPx(:,pindx.lbeta)];
            
            RHS(:,kk) = -tmp;
            kk = kk + 1;
        end
        Pxx = mfactor(FFp, RHS);
        toc
    end
end

