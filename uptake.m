function [G,Gx,Gxx] = uptake(par, parm)
% unpack the parameters to be optimized
on = true; off = false;
alpha = parm.alpha;
beta  = parm.beta;

dLdbeta = parm.dLdbeta;
% d2Ldbetadbeta = parm.d2Ldbetadbeta;
iwet = parm.iwet;
nwet = parm.nwet;
% N:P of uptake operator
po4obs = parm.po4obs(iwet);
c2p    = parm.c2p;
% N uptake operator
L = diag(parm.L);  
dLdbeta = diag(parm.dLdbeta); 
% d2Ldbetadbeta = diag(parm.d2Ldbetadbeta);

DIP = parm.DIP;
G   = alpha.*c2p.*L.*DIP; 
Gp  = alpha.*c2p.*L;   

% gradient of uptake operator
nx  = parm.nx;
Gpx = zeros(nwet,nx);
Gpx(:,par.pindx.lalpha) = alpha.*c2p.*L;              % dGdlog_alpha
Gpx(:,par.pindx.lbeta)  = beta*alpha.*c2p.*dLdbeta;   % dGdlog_beta

% Gradient
% grad DIP
if nargout >1
    DIPx = zeros(nwet,nx);
    np = 0; % count the number of tunable parameters
    if (par.opt_sigma == on)
        np = np + 1;
        isigma = par.pindx.lsigma;
        DIPx(:,isigma) = parm.Px(1:nwet,isigma);
    end
    
    if (par.opt_kappa_dp == on)
        np = np + 1;
        ikappa_dp = par.pindx.lkappa_dp;
        DIPx(:,ikappa_dp) = parm.Px(1:nwet,ikappa_dp);
    end
    
    if (par.opt_slopep == on)
        np = np + 1;
        islopep = par.pindx.lslopep;
        DIPx(:,islopep) = parm.Px(1:nwet,islopep);
    end
    
    if (par.opt_interpp == on)
        np = np + 1;
        iinterpp = par.pindx.linterpp;
        DIPx(:,iinterpp) = parm.Px(1:nwet,iinterpp);
        
    end
    
    if (par.opt_alpha == on)
        np = np + 1;
        ialpha = par.pindx.lalpha; 
        DIPx(:,ialpha) = parm.Px(1:nwet,ialpha);
    end
    
    if (par.opt_beta == on)
        np = np + 1;
        ibeta = par.pindx.lbeta; 
        DIPx(:,ibeta) = parm.Px(1:nwet,ibeta);
    end
    
    Gx = zeros(nwet,nx);
    Gx = d0(Gp)*DIPx+d0(DIP)*Gpx;
end

%% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
if (nargout >2)
    dSi2Cdaa = parm.dSi2Cdaa;
    dSi2Cdbb = parm.dSi2Cdbb;
    DIPxx = parm.Pxx(1:nwet,:);
    % Compute the hessian of the solution wrt the parameters
    kk = 1;
    DIPxx = zeros(nwet,n);
    % sigma sigma
    if (par.opt_sigma == on)
        kk = kk + 1;
    end

    % sigma kappa_dp
    if (par.opt_sigma == on & par.opt_kappa_dp == on)
        kk = kk + 1;
    end

    % sigma slopep
    if (par.opt_sigma == on & par.opt_slopep == on)
        kk = kk + 1;
    end

    % sigma interpp
    if (par.opt_sigma == on & par.opt_interpp == on)
        kk = kk + 1;
    end

    % sigma alpha
    if (par.opt_sigma == on & par.opt_alpha == on)
        kk = kk + 1;
    end

    % sigma beta
    if (par.opt_sigma == on & par.opt_beta == on)
        kk = kk + 1;
    end

    % kappa_dp kappa_dp
    if (par.opt_kappa_dp == on)
        kk = kk + 1;
    end

    % kappa_dp slopep
    if (par.opt_kappa_dp == on & par.opt_slopep == on)
        kk = kk + 1;
    end

    % kappa_dp interpp
    if (par.opt_kappa_dp == on & par.opt_interpp == on)
        kk = kk + 1;
    end

    % kappa_dp alpha
    if (par.opt_kappa_dp == on & par.opt_alpha == on)
        kk = kk + 1;
    end

    % kappa_dp beta
    if (par.opt_kappa_dp == on & par.opt_beta == on)
        kk = kk + 1;
    end

    % slopep slopep
    if (par.opt_slopep == on)
        kk = kk + 1;
    end
    
    % slopep interpp
    if (par.opt_slopep == on & par.opt_interpp == on)
        kk = kk + 1;
    end

    % slopep alpha
    if (par.opt_slopep == on & par.opt_alpha == on)
        kk = kk + 1;
    end

    % slopep beta
    if (par.opt_slopep == on & par.opt_beta == on)
        kk = kk + 1;
    end
    
    % interpp interpp
    if (par.opt_interpp == on)
        kk = kk + 1;
    end    

    % interpp alpha
    if (par.opt_interpp == on & par.opt_alpha == on)
        kk = kk + 1;
    end

    % interpp beta
    if (par.opt_interpp == on & par.opt_beta == on)
        kk = kk + 1;
    end
    
    % alpha alpha
    if (par.opt_alpha == on)
        d2Gpdlog_alpha2 = alpha*c2p.*L;
        Gxx(:,kk) = d2Gpdlog_alpha2.*DIP + ...
            2*Gpx(:, par.pindx.lalpha).*DIPx(:, par.pindx.lalpha) + ...
            Gp.*DIPxx(:, kk);
        kk = kk + 1;
    end

    % alpha beta
    if (par.opt_alpha == on & par.opt_beta == on)
        d2Gpdlog_alphadlogbeta = alpha*beta*dLdbeta.*c2p;
        Gxx(:,kk) = d2Gpdlog_alphadlogbeta .* DIP + ...
            Gpx(:,par.pindx.lalpha).*DIPx(par.pindx.lbeta) + ...
            Gpx(:,par.pindx.lbeta) .*DIPx(par.pindx.lalpha) + ...
            Gp.*DIPxx(:,kk);
        kk = kk + 1;
    end
    
    % beta beta
    if (par.opt_beta == on)
        d2Gpdlogbeta2 = beta*alpha*c2p.*dLdbeta+...
            beta^2*alpha*c2p.*d2Ldbetadbeta;
        Gxx(:,kk) = d2Gpdlogbeta2 .* DIP + ...
            2*Gpx(:, par.pindx.lbeta).*DIPx(:, par.pindx.lbeta) + ...
            Gp.*DIPxx(:, kk);
        kk = kk + 1;
    end
    % --------------------------------------------------------------
    % Si parameters
    % sigma at
    if (par.opt_sigma == on & par.opt_at == on)
        kk = kk + 1;
    end
    % sigma bt
    if (par.opt_sigma == on & par.opt_bt == on)
        kk = kk + 1;
    end
    % sigma aa
    if (par.opt_sigma == on & par.opt_aa == on)
        kk = kk + 1;
    end
    % sigma bb
    if (par.opt_sigma == on & par.opt_bb == on)
        kk = kk + 1;
    end
    % kappa_dp at
    if (par.opt_kappa_dp == on & par.opt_at == on)
        kk = kk + 1;
    end
    % kappa_dp bt
    if (par.opt_kappa_dp == on & par.opt_bt == on)
        kk = kk + 1;
    end
    % kappa_dp aa
    if (par.opt_kappa_dp == on & par.opt_aa == on)
        kk = kk + 1;
    end
    % kappa_dp bb
    if (par.opt_kappa_dp == on & par.opt_bb == on)

    end
    % slope at
    if (par.opt_slope == on & par.opt_at == on)
        kk = kk + 1;
    end
    % slope bt
    if (par.opt_slope == on & par.opt_bt == on)
        kk = kk + 1;
    end
    % slope aa
    if (par.opt_slope == on & par.opt_aa == on)
        kk = kk + 1;
    end
    % slope bb
    if (par.opt_slope == on & par.opt_bb == on)
        kk = kk + 1;
    end
    % interpp at
    if (par.opt_interpp == on & par.opt_at == on)
        kk = kk + 1;
    end
    % interpp bt
    if (par.opt_interpp == on & par.opt_bt == on)
        kk = kk + 1;
    end
    % interpp aa
    if (par.opt_interpp == on & par.opt_aa == on)
        kk = kk + 1;
    end
    % interpp bb
    if (par.opt_interpp == on & par.opt_bb == on)
        kk = kk + 1;
    end

    % alpha at
    if (par.opt_alpha == on & par.opt_at == on)
        kk = kk + 1;
    end
    % alpha bt
    if (par.opt_alpha == on & par.opt_bt == on)
        kk = kk + 1;
    end
    % alpha aa
    if (par.opt_alpha == on & par.opt_aa == on)
        Gxx(:,kk) = ...
            Gpx(:,par.pindx.lalpha).*DIP.*dSi2Cdaa + ...
            Gp.*DIPx(:,par.pindx.lalpha).*dSi2Cdaa;
        kk = kk + 1;
    end
    % alpha bb
    if (par.opt_alpha == on & par.opt_bb == on)
        Gxx(:,kk) = ...
            Gpx(:,par.pindx.lalpha).*DIP.*dSi2Cdbb + ...
            Gp.*DIPx(:,par.pindx.lalpha).*dSi2Cdbb;
        kk = kk + 1;
    end
    % beta at
    if (par.opt_beta == on & par.opt_at == on)
        kk = kk + 1;
    end
    % beta bt
    if (par.opt_beta == on & par.opt_bt == on)
        kk = kk + 1;
    end
    % beta aa
    if (par.opt_beta == on & par.opt_aa == on)
        Gxx(:,kk) = ...
            Gpx(:,par.pindx.lbeta).*DIP.*dSi2Cdaa + ...
            Gp.*DIPx(:,par.pindx.lbeta).*dSi2Cdaa;
        kk = kk + 1;
    end
    % beta bb
    if (par.opt_beta == on & par.opt_bb == on)
        Gxx(:,kk) = ...
            Gpx(:,par.pindx.lbeta).*DIP.*dSi2Cdbb + ...
            Gp.*DIPx(:,par.pindx.lbeta).*dSi2Cdbb;
        kk = kk + 1;
    end
    % at at
    if (par.opt_at == on)
        kk = kk + 1;
    end
    % at bt
    if (par.opt_at == on & par.opt_bt == on)
        kk = kk + 1;
    end
    % at aa
    if (par.opt_at == on & par.opt_aa == on)
        kk = kk + 1;
    end
    % at bb
    if (par.opt_at == on & par.opt_bb == on)
        kk = kk + 1;
    end
    % bt bt
    if (par.opt_bt == on)
        kk = kk + 1;
    end
    % bt aa
    if (par.opt_bt == on & par.opt_aa == on)
        kk = kk + 1;
    end
    % bt bb
    if (par.opt_bt == on & par.opt_bb == on)
        kk = kk + 1;
    end
    % aa aa
    if (par.opt_aa == on & par.opt_aa == on)
        kk = kk + 1;
    end
    % aa bb
    if (par.opt_aa == on & par.opt_bb == on)
        kk = kk + 1;
    end
    % bb bb
    if (par.opt_bb == on & par.opt_bb == on)
        Gxx(:,kk) = Gp.*DIP*.bb;
        kk = kk + 1;
    end
    
end

    % Hessian
    % map P-model parameter order to N-model parameter order

% n = (nx+1)*nx/2;
% pos = tril(ones(nx),0);
% pos(pos==1) = 1:n;
% pos = pos';

% grad grad DIP
% DIPxx = zeros(nwet,n);
% nn = 1;
% for ji = 1:np
    % for kk = ji:np
        % map(nn) = pos(ji,kk);
        % nn = nn + 1;
    % end
% end

% DIPxx(:,map) = parm.Pxx(1:nwet,:);

% Gpxx = zeros(nwet,n);
% Gpxx(:,pos(ialpha,ialpha)) = alpha*c2p.*L; % alpha alpha
% Gpxx(:,pos(ialpha,ibeta)) = beta*alpha*c2p.*dLdbeta; % alpha beta
% Gpxx(:,pos(3,9)) = c2p_l*alpha*dc2pdc2p_l.*L;  % alpha c2p_l 
% Gpxx(:,pos(3,10)) = S*alpha*dc2pdS.*L;  % alpha S 
% Gpxx(:,pos(ibeta,ibeta)) = beta*alpha*c2p.*dLdbeta+...
    % beta^2*alpha*c2p.*d2Ldbetadbeta;   % beta beta 
% Gpxx(:,pos(4,9)) = c2p_l*beta*alpha*dc2pdc2p_l.*dLdbeta; % beta c2p_l
% Gpxx(:,pos(4,10)) = S*beta*alpha*dc2pdS.*dLdbeta;         % beta S
% Gpxx(:,pos(9,9)) = c2p_l*alpha*dc2pdc2p_l.*L;           % c2p_l c2p_l
% Gpxx(:,pos(9,10)) = 0*G;                                 % c2p_l S    
% Gpxx(:,pos(10,10)) = S*alpha*dc2pdS.*L;                   % S S


% Gxx = d0(DIP)*Gpxx+d0(Gp)*DIPxx;
% for is = 1:nx
    % r = [pos(is,is):pos(is,end)];
    % Gxx(:,r) = Gxx(:,r)+...
        % d0(Gpx(:,is))*DIPx(:,is:end)+...
        % d0(DIPx(:,is))*Gpx(:,is:end);
% end





