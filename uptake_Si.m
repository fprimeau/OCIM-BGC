function [G,Gx,Gxx,Gp,par] = uptake_Si(par)
% unpack the parameters to be optimized
on = true; off = false;
iwet  = par.iwet;
nwet  = par.nwet;
alpha = par.alpha;
beta  = par.beta;

dLdbeta = diag(par.dLdbeta); 
d2Ldbetadbeta = diag(par.d2Ldbetadbeta);
% N:P of uptake operator
po4obs = par.po4obs(iwet);
c2p    = par.c2p;
Si2C   = par.Si2C;
% N uptake operator
L = diag(par.L);

DIP = par.DIP;
G   = d0(alpha.*c2p.*L.*DIP);
Gp  = d0(alpha.*c2p.*L);

% gradient of uptake operator
nx  = par.nx;
Gpx = zeros(nwet,nx);
Gpx(:,par.pindx.lalpha) = alpha.*c2p.*L;              % dGdlog_alpha
Gpx(:,par.pindx.lbeta)  = beta*alpha.*c2p.*dLdbeta;   % dGdlog_beta

% Gradient
% grad DIP
if (nargout >1)
    DIPx = zeros(nwet,nx);
    np = 0; % count the number of tunable parameters
    if (par.opt_sigma == on)
        isigma = par.pindx.lsigma;
        DIPx(:,isigma) = par.Px(1:nwet,isigma);
    end

    if (par.opt_kappa_dp == on)
        ikappa_dp = par.pindx.lkappa_dp;
        DIPx(:,ikappa_dp) = par.Px(1:nwet,ikappa_dp);
    end

    if (par.opt_slopep == on)
        islopep = par.pindx.slopep;
        DIPx(:,islopep) = par.Px(1:nwet,islopep);
    end

    if (par.opt_interpp == on)
        iinterpp = par.pindx.linterpp;
        DIPx(:,iinterpp) = par.Px(1:nwet,iinterpp);
    end

    if (par.opt_alpha == on)
        ialpha = par.pindx.lalpha;
        DIPx(:,ialpha) = par.Px(1:nwet,ialpha);
    end

    if (par.opt_beta == on)
        ibeta = par.pindx.lbeta;
        DIPx(:,ibeta) = par.Px(1:nwet,ibeta);
    end

    Gx = zeros(nwet,nx);
    Gx = Gp*DIPx+d0(DIP)*Gpx;

    if (par.Simodel)
        dSi2Cdaa = par.dSi2Cdaa;
        dSi2Cdbb = par.dSi2Cdbb;

        if par.opt_aa
            iaa = par.pindx.aa;
            Gx(:,iaa) = sparse(nwet,1);
        end

        if par.opt_bb
            ibb = par.pindx.lbb;
            Gx(:,ibb) = sparse(nwet,1);
        end 
    end
end

%% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
if (nargout >2)
    npx = par.npx;
    dGpdalpha = Gpx(:,par.pindx.lalpha);
    dGpdbeta = Gpx(:,par.pindx.lbeta);
    par.dGpdalpha = dGpdalpha;
    par.dGpdbeta  = dGpdbeta;
    DIPxx = par.Pxx(1:nwet,:);
    Gxx = sparse(nwet,nchoosek(npx,2)+npx);
    % Compute the hessian of the solution wrt the parameters
    kk = 1;
    % sigma sigma
    if (par.opt_sigma == on)
        Gxx(:,kk) = Gp*DIPxx(:, kk);
        kk = kk + 1;
    end
    % sigma kappa_dp
    if (par.opt_sigma == on & par.opt_kappa_dp == on)
        Gxx(:,kk) = Gp*DIPxx(:, kk);
        kk = kk + 1;
    end
    % sigma slopep
    if (par.opt_sigma == on & par.opt_slopep == on)
        Gxx(:,kk) = Gp*DIPxx(:, kk);
        kk = kk + 1;
    end
    % sigma interpp
    if (par.opt_sigma == on & par.opt_interpp == on)
        Gxx(:,kk) = Gp*DIPxx(:, kk);
        kk = kk + 1;
    end
    % sigma alpha
    if (par.opt_sigma == on & par.opt_alpha == on)
        Gxx(:,kk) = d0(dGpdalpha)*DIPx(:, par.pindx.lsigma) + ...
            Gp*DIPxx(:, kk);
        kk = kk + 1;
    end
    % sigma beta
    if (par.opt_sigma == on & par.opt_beta == on)
        Gxx(:,kk) = d0(dGpdbeta)*DIPx(:, par.pindx.lsigma) + ...
            Gp*DIPxx(:, kk);
        kk = kk + 1;
    end
    % kappa_dp kappa_dp
    if (par.opt_kappa_dp == on)
        Gxx(:,kk) = Gp*DIPxx(:, kk);
        kk = kk + 1;
    end
    % kappa_dp slopep
    if (par.opt_kappa_dp == on & par.opt_slopep == on)
        Gxx(:,kk) = Gp*DIPxx(:, kk);
        kk = kk + 1;
    end
    % kappa_dp interpp
    if (par.opt_kappa_dp == on & par.opt_interpp == on)
        Gxx(:,kk) = Gp*DIPxx(:, kk);
        kk = kk + 1;
    end
    % kappa_dp alpha
    if (par.opt_kappa_dp == on & par.opt_alpha == on)
        Gxx(:,kk) = d0(dGpdalpha)*DIPx(:, par.pindx.lkappa_dp) + ...
            Gp*DIPxx(:, kk);
        kk = kk + 1;
    end
    % kappa_dp beta
    if (par.opt_kappa_dp == on & par.opt_beta == on)
        Gxx(:,kk) = d0(dGpdbeta)*DIPx(:, par.pindx.lkappa_dp) + ...
            Gp*DIPxx(:, kk);
        kk = kk + 1;
    end
    % slopep slopep
    if (par.opt_slopep == on)
        Gxx(:,kk) = Gp*DIPxx(:, kk);
        kk = kk + 1;
    end
    % slopep interpp
    if (par.opt_slopep == on & par.opt_interpp == on)
        Gxx(:,kk) = Gp*DIPxx(:, kk);
        kk = kk + 1;
    end
    % slopep alpha
    if (par.opt_slopep == on & par.opt_alpha == on)
        Gxx(:, kk) = d0(dGpdalpha)*DIPx(:, par.pindx.slopep) + ...
            Gp*DIPxx(:, kk);
        kk = kk + 1;
    end
    % slopep beta
    if (par.opt_slopep == on & par.opt_beta == on)
        Gxx(:, kk) = d0(dGpdbeta)*DIPx(:, par.pindx.slopep) + ...
            Gp*DIPxx(:, kk);
        kk = kk + 1;
    end
    % interpp interpp
    if (par.opt_interpp == on)
        Gxx(:, kk) = Gp*DIPxx(:, kk);
        kk = kk + 1;
    end
    % interpp alpha
    if (par.opt_interpp == on & par.opt_alpha == on)
        Gxx(:, kk) = d0(dGpdalpha)*DIPx(:, par.pindx.linterpp) + ...
            Gp*DIPxx(:, kk);
        kk = kk + 1;
    end
    % interpp beta
    if (par.opt_interpp == on & par.opt_beta == on)
        Gxx(:, kk) = d0(dGpdbeta)*DIPx(:, par.pindx.linterpp) + ...
            Gp*DIPxx(:,kk);
        kk = kk + 1;
    end
    % alpha alpha
    if (par.opt_alpha == on)
        d2Gpdlog_alpha2 = alpha*c2p.*L;
        Gxx(:, kk) = d0(DIP)*d2Gpdlog_alpha2 + ...
            2*Gpx(:, par.pindx.lalpha).*DIPx(:, par.pindx.lalpha) + ...
            Gp*DIPxx(:, kk);
        kk = kk + 1;
    end
    % alpha beta
    if (par.opt_alpha == on & par.opt_beta == on)
        d2Gpdlog_alphadlogbeta = alpha*beta*dLdbeta.*c2p;
        Gxx(:, kk) = d0(DIP)*d2Gpdlog_alphadlogbeta + ...
            Gpx(:, par.pindx.lalpha).*DIPx(:, par.pindx.lbeta) + ...
            Gpx(:, par.pindx.lbeta) .*DIPx(:, par.pindx.lalpha) + ...
            Gp*DIPxx(:, kk);
        kk = kk + 1;
    end
    % beta beta
    if (par.opt_beta == on)
        d2Gpdlogbeta2 = beta*alpha*c2p.*dLdbeta+...
            beta^2*alpha*c2p.*d2Ldbetadbeta;
        Gxx(:, kk) = d0(DIP)*d2Gpdlogbeta2 + ...
            2*Gpx(:, par.pindx.lbeta).*DIPx(:, par.pindx.lbeta) + ...
            Gp*DIPxx(:, kk);
        kk = kk + 1;
    end
    % --------------------------------------------------------------
    % Si parameters
    if par.Simodel == on
        % ---------------------------------------------------------
        nsx = par.nsx;
        at = par.at;
        bt  = par.bt;
        aa  = par.aa;
        bb  = par.bb;
        bsi = par.bsi;
        kappa_gs  = par.kappa_gs;
        % number of parameter combinations;
        if (nsx >=2)
            ncs = nchoosek(nsx,2)+nsx + npx*nsx;
            Gxx = cat(2,Gxx, sparse(nwet,ncs));
        end
        %
        % ---------------------------------------------------------
        %
        % sigma bsi
        if (par.opt_sigma == on & par.opt_bsi == on)
            kk = kk + 1;
        end
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
            Gxx(:, kk) = Gp*DIPx(:,par.pindx.lsigma).*dSi2Cdaa;
            kk = kk + 1;
        end
        % sigma bb
        if (par.opt_sigma == on & par.opt_bb == on)
            Gxx(:, kk) = bb*Gp*DIPx(:,par.pindx.lsigma).*dSi2Cdbb;
            kk = kk + 1;
        end
        % kappa_dp bsi
        if (par.opt_kappa_dp == on & par.opt_bsi == on)
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
            Gxx(:, kk) = Gp*DIPx(:,par.pindx.lkappa_dp).*dSi2Cdaa;
            kk = kk + 1;
        end
        % kappa_dp bb
        if (par.opt_kappa_dp == on & par.opt_bb == on)
            Gxx(:, kk) = bb*Gp*DIPx(:,par.pindx.lkappa_dp).*dSi2Cdbb;
            kk = kk + 1;
        end
        % slope bsi
        if (par.opt_slopep == on & par.opt_bsi == on)
            kk = kk + 1;
        end
        % slope at
        if (par.opt_slopep == on & par.opt_at == on)
            kk = kk + 1;
        end
        % slope bt
        if (par.opt_slopep == on & par.opt_bt == on)
            kk = kk + 1;
        end
        % slope aa
        if (par.opt_slopep == on & par.opt_aa == on)
            Gxx(:, kk) = Gp*DIPx(:,par.pindx.slopep).*dSi2Cdaa;
            kk = kk + 1;
        end
        % slope bb
        if (par.opt_slopep == on & par.opt_bb == on)
            Gxx(:, kk) = bb*Gp*DIPx(:,par.pindx.slopep).*dSi2Cdbb;
            kk = kk + 1;
        end
        % interpp bsi
        if (par.opt_interpp == on & par.opt_bsi == on)
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
            Gxx(:, kk) = Gp*DIPx(:,par.pindx.linterpp).*dSi2Cdaa;
            kk = kk + 1;
        end
        % interpp bb
        if (par.opt_interpp == on & par.opt_bb == on)
            Gxx(:, kk) = bb*Gp*DIPx(:,par.pindx.linterpp).*dSi2Cdbb;
            kk = kk + 1;
        end
        % alpha bsi
        if (par.opt_alpha == on & par.opt_bsi == on)
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
            Gxx(:, kk) = ...
                d0(Gpx(:,par.pindx.lalpha))*DIP.*dSi2Cdaa + ...
                Gp*DIPx(:,par.pindx.lalpha).*dSi2Cdaa;
            kk = kk + 1;
        end
        % alpha bb
        if (par.opt_alpha == on & par.opt_bb == on)
            Gxx(:, kk) = bb * ...
                (d0(Gpx(:,par.pindx.lalpha))*DIP.*dSi2Cdbb + ...
                 Gp*DIPx(:,par.pindx.lalpha).*dSi2Cdbb);
            kk = kk + 1;
        end
        % beta bsi
        if (par.opt_beta == on & par.opt_bsi == on)
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
            Gxx(:, kk) = ...
                d0(Gpx(:,par.pindx.lbeta))*DIP.*dSi2Cdaa + ...
                Gp*DIPx(:,par.pindx.lbeta).*dSi2Cdaa;
            kk = kk + 1;
        end
        % beta bb
        if (par.opt_beta == on & par.opt_bb == on)
            Gxx(:, kk) = bb * ...
                (d0(Gpx(:,par.pindx.lbeta))*DIP.*dSi2Cdbb + ...
                 Gp*DIPx(:,par.pindx.lbeta).*dSi2Cdbb);
            kk = kk + 1;
        end
        % bsi bsi
        if (par.opt_bsi)
            kk = kk + 1;
        end
        % bsi at
        if (par.opt_bsi & par.opt_at)
            kk = kk + 1;
        end
        % bsi bt
        if (par.opt_bsi & par.opt_bt)
            kk = kk + 1;
        end
        % bsi aa
        if (par.opt_bsi & par.opt_aa)
            kk = kk + 1;
        end
        % bsi bb
        if (par.opt_bsi & par.opt_bb)
            kk = kk + 1;
        end
        % at at
        if (par.opt_at)
            kk = kk + 1;
        end
        % at bt
        if (par.opt_at & par.opt_bt)
            kk = kk + 1;
        end
        % at aa
        if (par.opt_at & par.opt_aa)
            kk = kk + 1;
        end
        % at bb
        if (par.opt_at & par.opt_bb)
            kk = kk + 1;
        end
        % bt bt
        if (par.opt_bt)
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
            Gxx(:, kk) = bb*G*dSi2Cdbb;
            kk = kk + 1;
        end 
    end
end

end

