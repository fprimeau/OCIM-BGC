function [G,Gx,Gxx] = production_O(par,parm)
% unpack the parameters to be optimized
alpha = parm.alpha;
beta  = parm.beta;

iwet = parm.iwet;
nwet = parm.nwet;
% of uptake operator
po4obs = parm.po4obs(iwet);
c2p = parm.c2p;
% P uptake operator
L = parm.L;  
dLdbeta = parm.dLdbeta; 


DIP = parm.DIP;
G   = d0(alpha*L*DIP); 
Gp  = alpha*L;   

pindx = par.pindx;
% Gradient
% grad DIP
if (nargout >1)
    % gradient of uptake operator
    nx  = parm.npx;
    Gpx = zeros(nwet,nx);

    DIPx = zeros(nwet,nx);
    np = 0; % count the number of tunable parameters
    if (par.opt_sigma)
        isigma = pindx.lsigma;
        DIPx(:,isigma) = parm.Px(1:nwet,isigma);
    end

    if (par.opt_kappa_dp)
        ikappa_dp = pindx.lkappa_dp;
        DIPx(:,ikappa_dp) = parm.Px(1:nwet,ikappa_dp);
    end

    if (par.opt_slopep)
        islopep = pindx.slopep;
        DIPx(:,islopep) = parm.Px(1:nwet,islopep);
    end

    if (par.opt_interpp)
        iinterpp = pindx.linterpp;
        DIPx(:,iinterpp) = parm.Px(1:nwet,iinterpp);
    end

    if (par.opt_alpha)
        ialpha = pindx.lalpha;
        Gpx(:,ialpha) = diag(alpha*L); % dGdlog_alpha
        DIPx(:,ialpha) = parm.Px(1:nwet,ialpha);
    end

    if (par.opt_beta)
        ibeta = pindx.lbeta;
        Gpx(:,ibeta) = diag(beta*alpha*dLdbeta); % dGdlog_beta
        DIPx(:,ibeta) = parm.Px(1:nwet,ibeta);
    end

    Gx = zeros(nwet,nx);
    Gx = Gp*DIPx+d0(DIP)*Gpx;
end

%% ------------------------------------------------
if (nargout >2)
    kk = 0;
    DIPxx = parm.Pxx(1:nwet,:);
    d2Ldbetadbeta = parm.d2Ldbetadbeta;
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
        d2Gpdlogbeta2 = beta*alpha*dLdbeta+...
            beta^2*alpha*diag(d2Ldbetadbeta);
        Gxx(:,kk) = Gp*DIPxx(:,kk) + ...
            d0(Gpx(:,pindx.lalpha))*DIPx(:,pindx.lbeta) + ...
            d2Gpdlogbeta2*DIP;
    end
end