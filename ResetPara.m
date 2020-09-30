function  x = ResetPara(x, par);
    global iter
    on = true; off = false;
    load(par.fxhat,'x0')
    %
    fb = 2 ;  
    % gradually decrease the constrains    
    fb = fb * 1.25^iter ;
    fs = 1./fb         ;
    fprintf('current fb is %3.3e \n', fb)
    % parameter values
    pindx = par.pindx  ;
    
    if (par.opt_sigma == on)
        isigma = pindx.lsigma      ;
        xnew   = exp(x(isigma))    ;
        xold   = exp(x0(isigma))   ;
        if (xnew > 1 | xnew < 0)
            x(isigma) = x0(isigma) ;
        end
    end

    if (par.opt_kP_T == on & par.opt_kdP == off)
        ikP_T = pindx.kP_T;
        kP_T1 = x(ikP_T)  ;
        kP_T0 = x0(ikP_T) ;
        %
        kP1  = kP_T1*par.Tz + par.kdP ;
        kP0  = kP_T0*par.Tz + par.kdP ;
        mkP1 = mean(kP1)  ;
        mkP0 = mean(kP0)  ; 
        if any(kP1 < 0) | mkP1 > fb*mkP0 | mkP1 < fs*mkP0
            x(ikP_T) = x0(ikP_T) ;
        end
    end

    if (par.opt_kP_T == on & par.opt_kdP == on)
        ikP_T = pindx.kP_T   ;
        kP_T1 = x(ikP_T)     ;
        kP_T0 = x0(ikP_T)    ;
        %
        ikdP = pindx.lkdP    ;
        kdP1 = exp(x(ikdP))  ;
        kdP0 = exp(x0(ikdP)) ;
        %
        kP1  = kP_T1*par.Tz + kdP1 ;
        kP0  = kP_T0*par.Tz + kdP0 ;
        mkP1 = mean(kP1) ;
        mkP0 = mean(kP0) ; 
        if any(kP1 < 0) | mkP1 > fb*mkP0 | mkP1 < fs*mkP0
            x(ikP_T) = x0(ikP_T) ;
            x(ikdP)  = x0(ikdP)  ;
        end
    end
    
    if (par.opt_kP_T == off & par.opt_kdP == on)
        ikdP = pindx.lkdP     ;
        xnew = exp(x(ikdP))   ;
        xold = exp(x0(ikdP))  ;
        if (xnew > fb*xold | xnew < fs*xold);
            x(ikdP) = x0(ikdP) ;
        end
    end

    if (par.opt_bP_T == on & par.opt_bP == off)
        ibP_T = pindx.bP_T ;
        bm1   = x(ibP_T)   ;
        bm0   = x0(ibP_T)  ;
        %
        bb1   = par.bP     ; 
        bP    = bm1*par.aveT + bb1 ;
        if min(bP(:)) < 0.3 | max(bP(:)) > 4
            x(ibP_T) = x0(ibP_T)   ;
        end
    end

    if (par.opt_bP_T == on & par.opt_bP == on)
        ibP_T = pindx.bP_T  ;
        bm1   = x(ibP_T)    ;
        bm0   = x0(ibP_T)   ;
        %
        ibP  = pindx.lbP    ;
        bb1  = exp(x(ibP))  ;
        bb0  = exp(x0(ibP)) ;
        %
        bP   = bm1*par.aveT + bb1 ;
        if min(bP(:)) < 0.3 | max(bP(:)) > 4
            x(ibP_T) = x0(ibP_T) ;
            x(ibP)   = x0(ibP);
        end
    end

    if (par.opt_bP_T == off & par.opt_bP == on)
        ibP  = pindx.lbP    ;
        xnew = exp(x(ibP))  ;
        xold = exp(x0(ibP)) ;
        if (xnew < 0.3 | xnew > 4) ;
            x(ibP) = x0(ibP);
        end
    end
    
    if (par.opt_alpha == on)
        ialpha = pindx.lalpha    ;
        xnew   = exp(x(ialpha))  ;
        xold   = exp(x0(ialpha)) ;
        if (xnew > fb*xold | xnew < fs*xold);
            x(ialpha) = x0(ialpha) ;
        end
    end

    if (par.opt_beta == on)
        ibeta = pindx.lbeta    ;
        xnew  = exp(x(ibeta))  ;
        xold  = exp(x0(ibeta)) ;
        if (xnew > fb*xold | xnew < fs*xold);
            x(ibeta) = x0(ibeta) ;
        end
    end

    if (par.Cmodel == on) 
        if (par.opt_bC_T == on & par.opt_bC == off)
            ibC_T = pindx.bC_T  ;
            bm1   = x(ibC_T)    ;
            bm0   = x0(ibC_T)   ;
            %
            bC   = bm1*par.aveT + par.bC ;
            if min(bC(:)) < 0.3 | max(bC(:)) > 4
                x(ibC_T) = x0(ibC_T) ;
            end
        end

        if (par.opt_bC_T == on  & par.opt_bC == on)
            ibC_T = pindx.bC_T  ;
            bm1   = x(ibC_T)    ;
            bm0   = x0(ibC_T)   ;
            %
            ibC  = pindx.lbC    ;
            bb1  = exp(x(ibC))  ;
            bb0  = exp(x0(ibC)) ;
            %
            bC   = bm1*par.aveT + bb1 ;
            if min(bC(:)) < 0.3 | max(bC(:)) > 4
                x(ibC_T) = x0(ibC_T) ;
                x(ibC)   = x0(ibC) ;
            end
        end
        
        if (par.opt_bC_T == off & par.opt_bC == on)
            ibC  = pindx.lbC    ;
            xnew = exp(x(ibC))  ;
            xold = exp(x0(ibC)) ;
            if (xnew < 0.3 | xnew > 4) ;
                x(ibC) = x0(ibC) ;
            end
        end

        if (par.opt_kPIC == on)
            ikPIC = pindx.lkPIC    ;
            xnew  = exp(x(ikPIC))  ;
            xold  = exp(x0(ikPIC)) ;
            if (xnew > fb*xold | xnew < fs*xold) ;
                x(ikPIC) = x0(ikPIC) ;
            end
        end
        
        if (par.opt_d == on)
            id   = pindx.ld    ;
            xnew = exp(x(id))  ;
            xold = exp(x0(id)) ;
            if (xnew > fb*xold | xnew < fs*xold) ;
                x(id) = x0(id) ;
            end
        end

        if (par.opt_kC_T == on & par.opt_kdC == off)
            ikC_T = pindx.kC_T ;
            kC_T1 = x(ikC_T)   ;
            kC_T0 = x0(ikC_T)  ;
            %
            kC1 = kC_T1*par.Tz + par.kdC ;
            kC0 = kC_T0*par.Tz + par.kdC ;
            mkC1 = mean(kC1) ;
            mkC0 = mean(kC0) ; 
            if any(kC1 < 0) | mkC1 > fb*mkC0 | mkC1 < fs*mkC0
                x(ikC_T) = x0(ikC_T);
            end
        end

        if (par.opt_kC_T == on & par.opt_kdC == on)
            ikC_T = pindx.kC_T ;
            kC_T1 = x(ikC_T)   ;
            kC_T0 = x0(ikC_T)  ;
            %
            ikdC  = pindx.lkdC    ;
            kdC1  = exp(x(ikdC))  ;
            kdC0  = exp(x0(ikdC)) ;
            %
            kC1 = kC_T1*par.Tz + kdC1 ;
            kC0 = kC_T0*par.Tz + kdC0 ;
            mkC1 = mean(kC1) ;
            mkC0 = mean(kC0) ; 
            if any(kC1 < 0) | mkC1 > fb*mkC0 | mkC1 < fs*mkC0
                x(ikC_T) = x0(ikC_T) ;
                x(ikdC)  = x0(ikdC) ;
            end
        end

        if (par.opt_kC_T == off & par.opt_kdC == on)
            ikdC = pindx.lkdC ;
            xnew = exp(x(ikdC))  ;
            xold = exp(x0(ikdC)) ;
            if (xnew > fb*xold | xnew < fs*xold) ;
                x(ikdC) = x0(ikdC) ;
            end
        end
        
        if (par.opt_R_Si == on & par.opt_rR == on)
            iR_Si = pindx.lR_Si ;
            irR   = pindx.lrR   ;
            
            par.R_Si = exp(x(iR_Si)) ;
            par.rR   = exp(x(irR))   ;
            vout  = mkPIC2P(par)     ;
            RR    = vout.RR          ;
            nRR   = mean(diag(RR))   ;

            par.R_Si = exp(x0(iR_Si)) ;
            par.rR   = exp(x0(irR))   ;
            vout  = mkPIC2P(par)      ;
            RR    = vout.RR           ;
            oRR   = mean(diag(RR))    ;
            
            if (nRR > fb*oRR | nRR < fs*oRR);
                x(iR_Si) = x0(iR_Si) ;
                x(irR)   = x0(irR)   ;
            end
        end
        
        if (par.opt_R_Si == on & par.opt_rR == off)
            iR_Si = pindx.lR_Si    ;
            
            par.R_Si = exp(x(iR_Si)) ;
            vout  = mkPIC2P(par)     ;
            RR    = vout.RR          ;
            nRR   = mean(diag(RR))   ;

            par.R_Si = exp(x0(iR_Si)) ;
            vout  = mkPIC2P(par)      ;
            RR    = vout.RR           ;
            oRR   = mean(diag(RR))    ;
            
            if (nRR > fb*oRR | nRR < fs*oRR);
                x(iR_Si) = x0(iR_Si) ;
            end
        end

        if (par.opt_rR == on & par.opt_R_Si == off)
            irR  = pindx.lrR    ;
            xnew = exp(x(irR))  ;
            xold = exp(x0(irR)) ;
            if (xnew > fb*xold  | xnew < fs*xold);
                x(irR) = x0(irR);
            end
        end
        
        if (par.opt_cc == on)
            icc  = pindx.lcc    ;
            xnew = exp(x(icc))  ;
            xold = exp(x0(icc)) ;
            if (xnew > fb*xold  | xnew < fs*xold);
                x(icc) = x0(icc);
            end
        end
        
        if (par.opt_dd == on)
            idd  = pindx.ldd    ;
            xnew = exp(x(idd))  ;
            xold = exp(x0(idd)) ;
            if (xnew > fb*xold  | xnew < fs*xold);
                x(idd) = x0(idd);
            end
        end
    end 
    % ------------------------------------------------------------
    if (par.Omodel == on)
        if (par.opt_O2C_T == on & par.opt_rO2C == off)
            iO2C_T = pindx.O2C_T ;
            O2C_T1 = x(iO2C_T)   ;
            O2C_T0 = x0(iO2C_T)  ;
            %
            O2C = O2C_T1*par.Tz + par.rO2C1 ;
            if (min(O2C) < 0)    ;
                x(iO2C_T) = x0(iO2C_T) ;
            end
        end

        if (par.opt_O2C_T == on & par.opt_rO2C == on)
            iO2C_T = pindx.O2C_T  ;
            O2C_T1 = x(iO2C_T)    ;
            O2C_T0 = x0(iO2C_T)   ;
            %
            irO2C = pindx.lrO2C   ;
            rO2C1 = exp(x(irO2C)) ;
            rO2C0 = exp(x0(irO2C));
            %
            O2C = O2C_T1*par.Tz*1e8 + rO2C1 ;
            if (min(O2C) < 0) 
                x(irO2C)  = x0(irO2C)  ;
                x(iO2C_T) = x0(iO2C_T) ;
            end
        end

        if (par.opt_O2C_T == off & par.opt_rO2C == on)
            irO2C = pindx.lrO2C    ;
            xnew  = exp(x(irO2C))  ;
            xold  = exp(x0(irO2C)) ;
            if (xnew > fb*xold | xnew < fs*xold);
                x(irO2C) = x0(irO2C);
            end
        end

        if (par.opt_O2P_T == on)
            iO2P_T = pindx.O2P_T        ;
            xnew   = abs(x(iO2P_T))     ;
            xold   = abs(x0(iO2P_T))    ;
            if (xnew > 50 | xnew < -50) ;
                x(iO2P_T) = x0(iO2P_T)  ;
            end
        end
        
        if (par.opt_rO2P == on)
            irO2P = pindx.lrO2P       ;
            xnew     = exp(x(irO2P))  ;
            xold     = exp(x0(irO2P)) ;
            if (xnew > fb*xold | xnew < fs*xold);
                x(irO2P) = x0(irO2P)  ;
            end
        end
    end
    % ------------------------------------------------------------
    % dsi
    if (par.Simodel==on)
        if (par.opt_dsi == on)
            idsi = pindx.ldsi    ;
            xnew = exp(x(idsi))  ;
            xold = exp(x0(idsi)) ;
            if (xnew > 5000 | xnew < 1000);
                x(idsi) = x0(idsi);
            end
        end
        % at
        if (par.opt_at == on)
            iat  = pindx.lat    ;
            xnew = exp(x(iat))  ;
            xold = exp(x0(iat)) ;
            if (xnew > fb*xold  | xnew < fs*xold);
                x(iat) = x0(iat);
            end
        end
        % bt
        if (par.opt_bt == on)
            ibt  = pindx.lbt     ;
            xnew = exp(x(ibt))   ;
            xold = exp(x0(ibt))  ;
            if (xnew > 15000     | xnew < 10000);
                x(ibt) = x0(ibt) ;
            end
        end
        % aa
        if (par.opt_aa == on)
            iaa  = pindx.aa   ;
            xnew = x(iaa)     ;
            xold = x0(iaa)    ;
            if (xnew > 50 | xnew < -50);
                x(iaa) = x0(iaa);
            end
        end
        % bb
        if (par.opt_bb == on)
            ibb  = pindx.lbb    ;
            xnew = exp(x(ibb))  ;
            xold = exp(x0(ibb)) ;
            if (xnew > fb*xold  | xnew < fs*xold);
                x(ibb) = x0(ibb);
            end
        end    
    end

end

