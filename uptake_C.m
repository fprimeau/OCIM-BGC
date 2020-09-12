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
        npx  = par.npx;
        Gpx  = zeros(nwet,npx);
        DIPx = par.Px(1:nwet,:);
        
        if (par.opt_alpha)
            ialpha = pindx.lalpha;
            Gpx(:,ialpha) = diag(alpha*L); % dGdlog_alpha
        end

        if (par.opt_beta)
            dLdbeta = par.dLdbeta;    
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
        
        % sigma sigma
        if (par.opt_sigma)
            kk = kk + 1;
            Gxx(:,kk) = Gp*DIPxx(:,kk);
        end

        % sigma kP_T
        if (par.opt_sigma & par.opt_kP_T)
            kk = kk + 1;
            Gxx(:,kk) = Gp*DIPxx(:,kk);
        end

        % sigma kdP
        if (par.opt_sigma & par.opt_kdP)
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
        %%%%
        % kP_T kP_T
        if (par.opt_kP_T)
            kk = kk + 1;
            Gxx(:,kk) = Gp*DIPxx(:,kk);
        end

        % kP_T kdP
        if (par.opt_kP_T & par.opt_kdP)
            kk = kk + 1;
            Gxx(:,kk) = Gp*DIPxx(:,kk);
        end

        % kP_T bP_T
        if (par.opt_kP_T & par.opt_bP_T)
            kk = kk + 1;
            Gxx(:,kk) = Gp*DIPxx(:,kk);
        end

        % kP_T bP
        if (par.opt_kP_T & par.opt_bP)
            kk = kk + 1;
            Gxx(:,kk) = Gp*DIPxx(:,kk);
        end

        % kP_T alpha
        if (par.opt_kP_T & par.opt_alpha)
            kk = kk + 1;
            Gxx(:,kk) = Gp*DIPxx(:, kk) + ...
                d0(Gpx(:,pindx.lalpha))*DIPx(:,pindx.kP_T);
        end

        % kP_T beta
        if (par.opt_kP_T & par.opt_beta)
            kk = kk + 1;
            Gxx(:,kk) = Gp*DIPxx(:,kk) + ...
                d0(Gpx(:,pindx.lbeta))*DIPx(:,pindx.kP_T);
        end
        
        % kdP kdP
        if (par.opt_kdP)
            kk = kk + 1;
            Gxx(:,kk) = Gp*DIPxx(:,kk);
        end

        % kdP bP_T
        if (par.opt_kdP & par.opt_bP_T)
            kk = kk + 1;
            Gxx(:,kk) = Gp*DIPxx(:,kk);
        end

        % kdP bP
        if (par.opt_kdP & par.opt_bP)
            kk = kk + 1;
            Gxx(:,kk) = Gp*DIPxx(:,kk);
        end

        % kdP alpha
        if (par.opt_kdP & par.opt_alpha)
            kk = kk + 1;
            Gxx(:,kk) = Gp*DIPxx(:,kk) + ...
                d0(Gpx(:,pindx.lalpha))*DIPx(:,pindx.lkdP);
        end

        % kdP beta
        if (par.opt_kdP & par.opt_beta)
            kk = kk + 1;
            Gxx(:,kk) = Gp*DIPxx(:,kk) + ...
                d0(Gpx(:,pindx.lbeta))*DIPx(:,pindx.lkdP);
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
            d2Ldbetadbeta = par.d2Ldbetadbeta;
            d2Gpdlogbeta2 = beta*alpha*dLdbeta + ...
                beta^2*alpha*d2Ldbetadbeta;
            Gxx(:,kk) = Gp*DIPxx(:,kk) + d2Gpdlogbeta2*DIP + ...
                2*d0(Gpx(:,pindx.lbeta))*DIPx(:,pindx.lbeta);
        end
    end
end

