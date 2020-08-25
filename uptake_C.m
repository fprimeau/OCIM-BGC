function [G,Gx,Gxx] = uptake_C(par)
on = true; off = false;
%
pindx = par.pindx;
% unpack the parameters to be optimized
alpha = par.alpha;
beta  = par.beta;

iwet = par.iwet;
nwet = par.nwet;
% of uptake operator
po4obs = par.po4obs(iwet);
% P uptake operator
L = par.L;  

DIP = par.DIP;
G   = d0(alpha*L*DIP); 
Gp  = alpha*L;   

% Gradient
% grad DIP
if par.optim == off
    Gx = [];
elseif (par.optim & nargout > 1)
    % gradient of uptake operator
    nx  = par.npx;
    Gpx = zeros(nwet,nx);
    DIPx = par.Px(1:nwet,:);
    dLdbeta = par.dLdbeta;    
    if (par.opt_alpha)
        ialpha = pindx.lalpha;
        Gpx(:,ialpha) = diag(alpha*L); % dGdlog_alpha
    end

    if (par.opt_beta)
        ibeta = pindx.lbeta;
        Gpx(:,ibeta) = diag(beta*alpha*dLdbeta); % dGdlog_beta
    end

    Gx = Gp*DIPx+d0(DIP)*Gpx;
end

%% ------------------------------------------------
if par.optim == off
    Gxx = [];
elseif (par.optim & nargout > 2)
    kk = 0;
    DIPxx = par.Pxx(1:nwet,:);
    d2Ldbetadbeta = par.d2Ldbetadbeta;
    % sigma sigma
    if (par.opt_sigma)
        kk = kk + 1;
        Gxx(:,kk) = Gp*DIPxx(:,kk);
    end

    % sigma kappa_dp
    if (par.opt_sigma & par.opt_kappa_dp)
        kk = kk + 1;
        Gxx(:,kk) = Gp*DIPxx(:,kk);
    end

    % sigma bP_T
    if (par.opt_sigma & par.opt_bP_T)
        kk = kk + 1;
        Gxx(:,kk) = Gp*DIPxx(:,kk);
    end

    % sigma bP
    if (par.opt_sigma & par.opt_bP)
        kk = kk + 1;
        Gxx(:,kk) = Gp*DIPxx(:,kk);
    end

    % sigma alpha
    if (par.opt_sigma & par.opt_alpha)
        kk = kk + 1;
        Gxx(:,kk) = Gp*DIPxx(:, kk) + ...
            d0(Gpx(:,pindx.lalpha))*DIPx(:,pindx.lsigma);
    end

    % sigma beta
    if (par.opt_sigma & par.opt_beta)
        kk = kk + 1;
        Gxx(:,kk) = Gp*DIPxx(:,kk) + ...
            d0(Gpx(:,pindx.lbeta))*DIPx(:,pindx.lsigma);
    end

    % kappa_dp kappa_dp
    if (par.opt_kappa_dp)
        kk = kk + 1;
        Gxx(:,kk) = Gp*DIPxx(:,kk);
    end

    % kappa_dp bP_T
    if (par.opt_kappa_dp & par.opt_bP_T)
        kk = kk + 1;
        Gxx(:,kk) = Gp*DIPxx(:,kk);
    end

    % kappa_dp bP
    if (par.opt_kappa_dp & par.opt_bP)
        kk = kk + 1;
        Gxx(:,kk) = Gp*DIPxx(:,kk);
    end

    % kappa_dp alpha
    if (par.opt_kappa_dp & par.opt_alpha)
        kk = kk + 1;
        Gxx(:,kk) = Gp*DIPxx(:,kk) + ...
            d0(Gpx(:,pindx.lalpha))*DIPx(:,pindx.lkappa_dp);
    end

    % kappa_dp beta
    if (par.opt_kappa_dp & par.opt_beta)
        kk = kk + 1;
        Gxx(:,kk) = Gp*DIPxx(:,kk) + ...
            d0(Gpx(:,pindx.lbeta))*DIPx(:,pindx.lkappa_dp);
    end

    % bP_T bP_T
    if (par.opt_bP_T)
        kk = kk + 1;
        Gxx(:,kk) = Gp*DIPxx(:,kk);
    end

    % bP_T bP
    if (par.opt_bP_T & par.opt_bP)
        kk = kk + 1;
        Gxx(:,kk) = Gp*DIPxx(:,kk);
    end

    % bP_T alpha
    if (par.opt_bP_T & par.opt_alpha)
        kk = kk + 1;
        Gxx(:,kk) = Gp*DIPxx(:,kk) + ...
            d0(Gpx(:,pindx.lalpha))*DIPx(:,pindx.bP_T);
    end

    % bP_T beta
    if (par.opt_bP_T & par.opt_beta)
        kk = kk + 1;
        Gxx(:,kk) = Gp*DIPxx(:,kk) + ...
            d0(Gpx(:,pindx.lbeta))*DIPx(:,pindx.bP_T);
    end

    % bP bP
    if (par.opt_bP)
        kk = kk + 1;
        Gxx(:,kk) = Gp*DIPxx(:,kk);
    end

    % bP alpha
    if (par.opt_bP & par.opt_alpha)
        kk = kk + 1;
        Gxx(:,kk) = Gp*DIPxx(:,kk) + ...
            d0(Gpx(:,pindx.lalpha))*DIPx(:,pindx.lbP);
    end

    % bP beta
    if (par.opt_bP & par.opt_beta)
        kk = kk + 1;
        Gxx(:,kk) = Gp*DIPxx(:,kk) + ...
            d0(Gpx(:,pindx.lbeta))*DIPx(:,pindx.lbP);
    end

    % alpha alpha
    if (par.opt_alpha)
        kk = kk + 1;        
        d2Gpdlog_alpha2 = alpha*L;
        Gxx(:,kk) = Gp*DIPxx(:,kk) + d2Gpdlog_alpha2*DIP + ...
            2*d0(Gpx(:,pindx.lalpha))*DIPx(:,pindx.lalpha);
    end

    % alpha beta
    if (par.opt_alpha & par.opt_beta)
        kk = kk + 1;
        d2Gpdlog_alphadlogbeta = alpha*beta*dLdbeta;
        Gxx(:,kk) = Gp*DIPxx(:,kk) + ...
            d0(Gpx(:,pindx.lalpha))*DIPx(:,pindx.lbeta) + ...
            d0(Gpx(:,pindx.lbeta))*DIPx(:,pindx.lalpha) + ...
            d2Gpdlog_alphadlogbeta*DIP;
    end

    % beta beta
    if (par.opt_beta)
        kk = kk + 1;
        d2Gpdlogbeta2 = beta*alpha*dLdbeta + ...
            beta^2*alpha*d2Ldbetadbeta;
        Gxx(:,kk) = Gp*DIPxx(:,kk) + d2Gpdlogbeta2*DIP + ...
            2*d0(Gpx(:,pindx.lbeta))*DIPx(:,pindx.lbeta);
    end
end