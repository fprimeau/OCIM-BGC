function [Par, Fe, Fex, Fexx] = eqFecycle(X, Par)
    global GFe
    on = true; off = false;
    iwet = par.iwet;
    nwet = par.nwet;

    % unpack the parameters to be optimized

    %
    options.iprint = 1     ;
    options.atol   = 1e-13 ;
    options.rtol   = 1e-13 ;
    
    fprintf('Solving Fe model ...\n') ;
    X0 = GFe;
    [Fe, ierr] = nsnew(X0, @(X) Fe_eqn(X, par), options) ;
    if (ierr ~=0)
        fprintf('O2model did not converge. \n') ;
        fprintf('Using the time stepping method for initial Fe ... \n')
        [~, ~, ~, ~, Fetstep] = Fe_eqn(X0, par) ;
        nsteps = size(Fetstep, 2) ;
        dsteps = 4 ;            % can adjust for choosing initial Fe results
        %
        for n = 1 : dsteps ;
            current_step = round(nstep*(n/dsteps)) ;
            fprintf('Used Fe(:,:,:,%d) for initial Fe field of nsnew ...\n', current_step)
            current_Fe = Fetstep(:, current_step) ;
            [Fe, ierr] = nsnew(current_Fe, @(X) Fe_eqn(X, par), options) ;
            %
            if (ierr == 0)
                fprintf('Time stepping method works for estimating inital values and the Fe model converged.\n')

                fprintf('reset the global variable for the next call eqFecycle. \n')
                GFe = real(Fe) + 1e-7*randn(par.nwet,1) ;
                [F, FD, Fex, Fexx] = Fe_eqn(Fe, par) ;
                break;
            end
        end
    else
        fprintf('reset the global variable for the next call eqFecycle. \n')
        GFe = real(Fe) + 1e-7*randn(par.nwet,1) ;
        [F, FD, Fex, Fexx] = Fe_eqn(Fe, par) ;
    end
end

function [F, FD, Fex, Fexx, Fetstep] = Fe_eqn(Fe, par)
    
    %unpack some useful stuff
    on = true; off = false;
    grd   = par.grd   ;
    M3d   = par.M3d   ;
    TRdiv = par.TRdiv ;
    iwet  = par.iwet  ;
    nwet  = par.nwet  ;
    dVt   = par.dVt   ;
    I     = par.I     ;
    Tz    = par.Tz    ;

    tFe = X(0*nwet+1:1*nwet) ; 

    % variables from C and O model

    % fixed parameters

    % tunable parameters

    % Fe function
    % -> Frantz 처럼 tFe로만 할 것이냐?
    % -> Roshan 처럼 FeOr과 FeOH로 나눌 것이냐?
    % -> Roshan의 FeOH는 advection, PFD, Burial 다함...
    tFe = frFe + LFe ;

    TRdiv*tFe = Source + Biological uptake + Sink(Scavenging) ; 
    
    
    
    
    
    % F = ...

    % KFe = FeL ./ freeFe * Lig
    KFe = 1e-10;
    FeL = KFe * (freeFe * Lig)

    eq1 = TRdiv*DFe %DFe

    eq2 = PFD*PFe %





    
    % factorizize Jacobian matrix
    if (nargout > 1)
        FD = mfactor(TRdiv + d0(dLdO) + KO2) ;
    end
    
    %% ----------------Gradient--------------------

    %% ----------------Hessian---------------------

    %% -----------Time stepping method-------------
    % Using a semi-implicit scheme for time stepping method.
    % A: Trapezoide rule, B: Euler Forward)
    %
    %  d
    % -- O2 +TRdiv*O2(t) = PO2 - LO2 .*R(t) + KO2*(O2sat - O2(t))  
    % dt
    %
    if nargout > 4
        % time-step size schedule
        %                  1 hr  1 day  1 month 1 year 4 years   10 years
        dt_size = spa./[ 365*24   365       12     1     0.25       0.1 ];  
        nsteps  =      [  100     100      100   200     200        300 ];  
    
        %
        Fetstep    = zeros(length(Fe), length(nsteps));
        j = 0; 
        t = 0;
        
        %
        for k = 1:length(dt_size) % loop over step size schedule
               % set the step size
             dt = dt_size(k);
             % set the number of time steps
             n = nsteps(k);
             % trapezoid rule
             A = I + (dt/2)*(TRdiv+KO2);
             B = I - (dt/2)*(TRdiv+KO2);
             FA = mfactor(A);
           
            fprintf('dt = %f nsteps = %i \n', dt, n);
             for i = 1:n   % loop over time steps
                 Fedt = mfactor(FA, (B*Fe + dt*(-LFe.*R + PFe +KFe*Fesat)));
                 t = t+dt;
                 j = j+1;
                 Fetstep(:,j)  = Fedt;
                 T(j) = t;
                 fprintf('.');
             end
             fprintf('\n');
        end
    end
end