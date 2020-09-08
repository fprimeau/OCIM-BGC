function [G,Gx,Gxx] = uptake_C(par)
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
dLdbeta = par.dLdbeta;

DIP = par.DIP;
G   = d0(alpha*L*DIP); 
Gp  = alpha*L;   

% Gradient
% grad DIP
if (nargout >1)
    % gradient of uptake operator
    nx  = par.npx;
    Gpx = zeros(nwet,nx);
    DIPx = par.Px(1:nwet,:);
    
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
if (nargout >2)
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

    % sigma slopep
    if (par.opt_sigma & par.opt_slopep)
        kk = kk + 1;
        Gxx(:,kk) = Gp*DIPxx(:,kk);
    end

    % sigma interpp
    if (par.opt_sigma & par.opt_interpp)
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

    % kappa_dp slopep
    if (par.opt_kappa_dp & par.opt_slopep)
        kk = kk + 1;
        Gxx(:,kk) = Gp*DIPxx(:,kk);
    end

    % kappa_dp interpp
    if (par.opt_kappa_dp & par.opt_interpp)
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

    % slopep slopep
    if (par.opt_slopep)
        kk = kk + 1;
        Gxx(:,kk) = Gp*DIPxx(:,kk);
    end

    % slopep interpp
    if (par.opt_slopep & par.opt_interpp)
        kk = kk + 1;
        Gxx(:,kk) = Gp*DIPxx(:,kk);
    end

    % slopep alpha
    if (par.opt_slopep & par.opt_alpha)
        kk = kk + 1;
        Gxx(:,kk) = Gp*DIPxx(:,kk) + ...
            d0(Gpx(:,pindx.lalpha))*DIPx(:,pindx.slopep);
    end

    % slopep beta
    if (par.opt_slopep & par.opt_beta)
        kk = kk + 1;
        Gxx(:,kk) = Gp*DIPxx(:,kk) + ...
            d0(Gpx(:,pindx.lbeta))*DIPx(:,pindx.slopep);
    end

    % interpp interpp
    if (par.opt_interpp)
        kk = kk + 1;
        Gxx(:,kk) = Gp*DIPxx(:,kk);
    end

    % interpp alpha
    if (par.opt_interpp & par.opt_alpha)
        kk = kk + 1;
        Gxx(:,kk) = Gp*DIPxx(:,kk) + ...
            d0(Gpx(:,pindx.lalpha))*DIPx(:,pindx.linterpp);
    end

    % interpp beta
    if (par.opt_interpp & par.opt_beta)
        kk = kk + 1;
        Gxx(:,kk) = Gp*DIPxx(:,kk) + ...
            d0(Gpx(:,pindx.lbeta))*DIPx(:,pindx.linterpp);
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
end

