function [G,Gx,Gxx,Gp,par] = uptake_Si(par)
% unpack the parameters to be optimized
    on = true; off = false;
    iwet  = par.iwet;
    nwet  = par.nwet;
    alpha = par.alpha;
    beta  = par.beta;
    pindx = par.pindx;
    dLdbeta = diag(par.dLdbeta); 
    d2Ldbetadbeta = diag(par.d2Ldbetadbeta);
    % N:P of uptake operator
    po4obs = par.po4obs(iwet);
    c2p    = par.c2p;
    % N uptake operator
    L = diag(par.L);

    DIP = par.DIP;
    G   = d0(alpha.*c2p.*L.*DIP);
    Gp  = d0(alpha.*c2p.*L);

    % gradient of uptake operator
    nx  = par.nx;
    Gpx = zeros(nwet,nx);
    Gpx(:,pindx.lalpha) = alpha.*c2p.*L;              % dGdlog_alpha
    Gpx(:,pindx.lbeta)  = beta*alpha.*c2p.*dLdbeta;   % dGdlog_beta

    % Gradient
    % grad DIP
    if (par.optim == off)
        Gx = [];
    elseif (par.optim & nargout > 1)
        DIPx = zeros(nwet,nx);
        np = 0; % count the number of tunable parameters
        if (par.opt_sigma == on)
            isigma = pindx.lsigma;
            DIPx(:,isigma) = par.Px(1:nwet,isigma);
        end

        if (par.opt_kdP == on)
            ikdP = pindx.lkdP;
            DIPx(:,ikdP) = par.Px(1:nwet,ikdP);
        end

        if (par.opt_bP_T == on)
            ibP_T = pindx.bP_T;
            DIPx(:,ibP_T) = par.Px(1:nwet,ibP_T);
        end

        if (par.opt_bP == on)
            ibP = pindx.lbP;
            DIPx(:,ibP) = par.Px(1:nwet,ibP);
        end

        if (par.opt_alpha == on)
            ialpha = pindx.lalpha;
            DIPx(:,ialpha) = par.Px(1:nwet,ialpha);
        end

        if (par.opt_beta == on)
            ibeta = pindx.lbeta;
            DIPx(:,ibeta) = par.Px(1:nwet,ibeta);
        end

        Gx = zeros(nwet,nx);
        Gx = Gp*DIPx+d0(DIP)*Gpx;

        if (par.Simodel)
            dSi2Cdaa = par.dSi2Cdaa;
            dSi2Cdbb = par.dSi2Cdbb;

            if par.opt_aa
                iaa = pindx.aa;
                Gx(:,iaa) = sparse(nwet,1);
            end

            if par.opt_bb
                ibb = pindx.lbb;
                Gx(:,ibb) = sparse(nwet,1);
            end 
        end
    end

    %% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    if (par.optim == off)
        Gxx = [];
    elseif (par.optim & nargout >2)
        npx = par.npx;
        dGpdalpha = Gpx(:,pindx.lalpha);
        dGpdbeta = Gpx(:,pindx.lbeta);
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
        % sigma kdP
        if (par.opt_sigma == on & par.opt_kdP == on)
            Gxx(:,kk) = Gp*DIPxx(:, kk);
            kk = kk + 1;
        end
        % sigma bP_T
        if (par.opt_sigma == on & par.opt_bP_T == on)
            Gxx(:,kk) = Gp*DIPxx(:, kk);
            kk = kk + 1;
        end
        % sigma bP
        if (par.opt_sigma == on & par.opt_bP == on)
            Gxx(:,kk) = Gp*DIPxx(:, kk);
            kk = kk + 1;
        end
        % sigma alpha
        if (par.opt_sigma == on & par.opt_alpha == on)
            Gxx(:,kk) = d0(dGpdalpha)*DIPx(:, pindx.lsigma) + ...
                Gp*DIPxx(:, kk);
            kk = kk + 1;
        end
        % sigma beta
        if (par.opt_sigma == on & par.opt_beta == on)
            Gxx(:,kk) = d0(dGpdbeta)*DIPx(:, pindx.lsigma) + ...
                Gp*DIPxx(:, kk);
            kk = kk + 1;
        end
        % kdP kdP
        if (par.opt_kdP == on)
            Gxx(:,kk) = Gp*DIPxx(:, kk);
            kk = kk + 1;
        end
        % kdP bP_T
        if (par.opt_kdP == on & par.opt_bP_T == on)
            Gxx(:,kk) = Gp*DIPxx(:, kk);
            kk = kk + 1;
        end
        % kdP bP
        if (par.opt_kdP == on & par.opt_bP == on)
            Gxx(:,kk) = Gp*DIPxx(:, kk);
            kk = kk + 1;
        end
        % kdP alpha
        if (par.opt_kdP == on & par.opt_alpha == on)
            Gxx(:,kk) = d0(dGpdalpha)*DIPx(:, pindx.lkdP) + ...
                Gp*DIPxx(:, kk);
            kk = kk + 1;
        end
        % kdP beta
        if (par.opt_kdP == on & par.opt_beta == on)
            Gxx(:,kk) = d0(dGpdbeta)*DIPx(:, pindx.lkdP) + ...
                Gp*DIPxx(:, kk);
            kk = kk + 1;
        end
        % bP_T bP_T
        if (par.opt_bP_T == on)
            Gxx(:,kk) = Gp*DIPxx(:, kk);
            kk = kk + 1;
        end
        % bP_T bP
        if (par.opt_bP_T == on & par.opt_bP == on)
            Gxx(:,kk) = Gp*DIPxx(:, kk);
            kk = kk + 1;
        end
        % bP_T alpha
        if (par.opt_bP_T == on & par.opt_alpha == on)
            Gxx(:, kk) = d0(dGpdalpha)*DIPx(:, pindx.bP_T) + ...
                Gp*DIPxx(:, kk);
            kk = kk + 1;
        end
        % bP_T beta
        if (par.opt_bP_T == on & par.opt_beta == on)
            Gxx(:, kk) = d0(dGpdbeta)*DIPx(:, pindx.bP_T) + ...
                Gp*DIPxx(:, kk);
            kk = kk + 1;
        end
        % bP bP
        if (par.opt_bP == on)
            Gxx(:, kk) = Gp*DIPxx(:, kk);
            kk = kk + 1;
        end
        % bP alpha
        if (par.opt_bP == on & par.opt_alpha == on)
            Gxx(:, kk) = d0(dGpdalpha)*DIPx(:, pindx.lbP) + ...
                Gp*DIPxx(:, kk);
            kk = kk + 1;
        end
        % bP beta
        if (par.opt_bP == on & par.opt_beta == on)
            Gxx(:, kk) = d0(dGpdbeta)*DIPx(:, pindx.lbP) + ...
                Gp*DIPxx(:,kk);
            kk = kk + 1;
        end
        % alpha alpha
        if (par.opt_alpha == on)
            d2Gpdlog_alpha2 = alpha*c2p.*L;
            Gxx(:, kk) = d0(DIP)*d2Gpdlog_alpha2 + ...
                2*Gpx(:, pindx.lalpha).*DIPx(:, pindx.lalpha) + ...
                Gp*DIPxx(:, kk);
            kk = kk + 1;
        end
        % alpha beta
        if (par.opt_alpha == on & par.opt_beta == on)
            d2Gpdlog_alphadlogbeta = alpha*beta*dLdbeta.*c2p;
            Gxx(:, kk) = d0(DIP)*d2Gpdlog_alphadlogbeta + ...
                Gpx(:, pindx.lalpha).*DIPx(:, pindx.lbeta) + ...
                Gpx(:, pindx.lbeta) .*DIPx(:, pindx.lalpha) + ...
                Gp*DIPxx(:, kk);
            kk = kk + 1;
        end
        % beta beta
        if (par.opt_beta == on)
            d2Gpdlogbeta2 = beta*alpha*c2p.*dLdbeta+...
                beta^2*alpha*c2p.*d2Ldbetadbeta;
            Gxx(:, kk) = d0(DIP)*d2Gpdlogbeta2 + ...
                2*Gpx(:, pindx.lbeta).*DIPx(:, pindx.lbeta) + ...
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
            dsi = par.dsi;
            kappa_g  = par.kappa_g;
            % number of parameter combinations;
            if (nsx >=2)
                ncs = nchoosek(nsx,2)+nsx + npx*nsx;
                Gxx = cat(2,Gxx, sparse(nwet,ncs));
            end
            %
            % ---------------------------------------------------------
            %
            % sigma dsi
            if (par.opt_sigma == on & par.opt_dsi == on)
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
                Gxx(:, kk) = Gp*DIPx(:,pindx.lsigma).*dSi2Cdaa;
                kk = kk + 1;
            end
            % sigma bb
            if (par.opt_sigma == on & par.opt_bb == on)
                Gxx(:, kk) = bb*Gp*DIPx(:,pindx.lsigma).*dSi2Cdbb;
                kk = kk + 1;
            end
            % kdP dsi
            if (par.opt_kdP == on & par.opt_dsi == on)
                kk = kk + 1;
            end
            % kdP at
            if (par.opt_kdP == on & par.opt_at == on)
                kk = kk + 1;
            end
            % kdP bt
            if (par.opt_kdP == on & par.opt_bt == on)
                kk = kk + 1;
            end
            % kdP aa
            if (par.opt_kdP == on & par.opt_aa == on)
                Gxx(:, kk) = Gp*DIPx(:,pindx.lkdP).*dSi2Cdaa;
                kk = kk + 1;
            end
            % kdP bb
            if (par.opt_kdP == on & par.opt_bb == on)
                Gxx(:, kk) = bb*Gp*DIPx(:,pindx.lkdP).*dSi2Cdbb;
                kk = kk + 1;
            end
            % slope dsi
            if (par.opt_bP_T == on & par.opt_dsi == on)
                kk = kk + 1;
            end
            % slope at
            if (par.opt_bP_T == on & par.opt_at == on)
                kk = kk + 1;
            end
            % slope bt
            if (par.opt_bP_T == on & par.opt_bt == on)
                kk = kk + 1;
            end
            % slope aa
            if (par.opt_bP_T == on & par.opt_aa == on)
                Gxx(:, kk) = Gp*DIPx(:,pindx.bP_T).*dSi2Cdaa;
                kk = kk + 1;
            end
            % slope bb
            if (par.opt_bP_T == on & par.opt_bb == on)
                Gxx(:, kk) = bb*Gp*DIPx(:,pindx.bP_T).*dSi2Cdbb;
                kk = kk + 1;
            end
            % bP dsi
            if (par.opt_bP == on & par.opt_dsi == on)
                kk = kk + 1;
            end
            % bP at
            if (par.opt_bP == on & par.opt_at == on)
                kk = kk + 1;
            end
            % bP bt
            if (par.opt_bP == on & par.opt_bt == on) 
                kk = kk + 1;
            end
            % bP aa
            if (par.opt_bP == on & par.opt_aa == on)
                Gxx(:, kk) = Gp*DIPx(:,pindx.lbP).*dSi2Cdaa;
                kk = kk + 1;
            end
            % bP bb
            if (par.opt_bP == on & par.opt_bb == on)
                Gxx(:, kk) = bb*Gp*DIPx(:,pindx.lbP).*dSi2Cdbb;
                kk = kk + 1;
            end
            % alpha dsi
            if (par.opt_alpha == on & par.opt_dsi == on)
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
                    d0(Gpx(:,pindx.lalpha))*DIP.*dSi2Cdaa + ...
                    Gp*DIPx(:,pindx.lalpha).*dSi2Cdaa;
                kk = kk + 1;
            end
            % alpha bb
            if (par.opt_alpha == on & par.opt_bb == on)
                Gxx(:, kk) = bb * ...
                    (d0(Gpx(:,pindx.lalpha))*DIP.*dSi2Cdbb + ...
                     Gp*DIPx(:,pindx.lalpha).*dSi2Cdbb);
                kk = kk + 1;
            end
            % beta dsi
            if (par.opt_beta == on & par.opt_dsi == on)
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
                    d0(Gpx(:,pindx.lbeta))*DIP.*dSi2Cdaa + ...
                    Gp*DIPx(:,pindx.lbeta).*dSi2Cdaa;
                kk = kk + 1;
            end
            % beta bb
            if (par.opt_beta == on & par.opt_bb == on)
                Gxx(:, kk) = bb * ...
                    (d0(Gpx(:,pindx.lbeta))*DIP.*dSi2Cdbb + ...
                     Gp*DIPx(:,pindx.lbeta).*dSi2Cdbb);
                kk = kk + 1;
            end
            % dsi dsi
            if (par.opt_dsi)
                kk = kk + 1;
            end
            % dsi at
            if (par.opt_dsi & par.opt_at)
                kk = kk + 1;
            end
            % dsi bt
            if (par.opt_dsi & par.opt_bt)
                kk = kk + 1;
            end
            % dsi aa
            if (par.opt_dsi & par.opt_aa)
                kk = kk + 1;
            end
            % dsi bb
            if (par.opt_dsi & par.opt_bb)
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

