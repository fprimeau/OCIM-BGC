function [PFdiv,Gout,Hout] = buildPFD(par,Ptype)
%[PFdiv] = buildPFD(M3d,grd,par,Ptype)
%[PFdiv,grad] = buildPFD(M3d,grd,par,Ptype)
%[PFdiv,grad,hess] = buildPFD(M3d,grd,par,Ptype)
%
% INPUT:
%   M3D: the wet=1, dry = 0 mask
%   grd: the grid definition structure
%   par: parameters and auxiliary variables
%   PTtype: sinking particle type
%           options: 'POC','POP','PIC','bSi'
%
% this code is used to build a Particle Flux Diverngence operator
%
%                      ________________________
%                     /         A             /|  POM sinking flux
%                 top/_______________________/ |       |
%                    |  |                    | |       | -w
%                    | /        V            | /       |
%                 bot|/______________________|/        V
%
    M3d = par.M3d;
    grd = par.grd;
    [DIV,IU,M,MSK,iwet] = mkOperators(M3d,grd);
    
    if ( strcmp(Ptype,'POP') | strcmp(Ptype,'POC') )
        % Define a temperature dependent powerlaw exponent
        T = par.aveT;
        if ( strcmp(Ptype,'POP') )
            % extract the parameters
            % use bb = 0 for a constant b value
            bb = par.bP;  % intercept
            bm = par.bP_T; % slope
        elseif ( strcmp(Ptype,'POC') )
            % extract the parameters
            % use bm = 0 for a constant b value
            bb = par.bC;   % intercept
            bm = par.bC_T; % slope
        end
        b = bb+bm*T;
        % particle sinking velocity at the top of the grid cells.
        % mimic a Martin curve flux attenuation profile (see Kriest and
        % Oschelies 2008 in Biogeosciences)
        r =  par.kappa_p; % dissolution rate
        a =  r./b;
        w = -a.*M;
        % 'upwind' sinking flux operator
        % defined at the top of the grid cell
        FLUX = d0(w(iwet))*IU;
    elseif ( strcmp(Ptype,'PIC') )
        r = 1./par.taup;
        d = par.d; % e-folding length scale for PIC flux attenuation
        w = -r*d*MSK;
        FLUX = d0(w(iwet))*IU;
    elseif ( strcmp(Ptype,'bSi') )
        % Define a temperature dependent dissolution rate
        % Sarmiento and Gruber text book
        T  =  par.modT + 273.15;
        at =  par.at; 
        bt =  par.bt;
        d  =  par.dsi;
        r  =  at*exp(-bt./T);
        w  =  MSK; 
        w(:,:,2:end-1)  = -r(:,:,2:end)*d;
        FLUX = d0(w(iwet))*IU;
    end
    % Particle flux divergence operator
    PFdiv = DIV*FLUX;

    if (nargout > 1) % compute the gradient w.r.t. the parameters
        if ( strcmp(Ptype,'POP') | strcmp(Ptype,'POC') )
            a_b  = -r./b.^2;  % derivative of a w.r.t. b
            a_r  = 1./b;      % derivative of a w.r.t. r
            b_bb = 1;         % derivative of b w.r.t. bb
            b_bm = T;         % derivative of b w.r.t. bm
            w_b  = -a_b.*M;   % derivative of w w.r.t. b
            w_r  = -a_r.*M;   % derivative of w w.r.t. r
            w_bm = w_b.*b_bm; % derivative of w w.r.t. bm
            w_bb = w_b.*b_bb; % derivative of w w.r.t. bb
            
            FLUX_r = d0(w_r(iwet))*IU;
            FLUX_b = d0(w_b(iwet))*IU;
            FLUX_bm = d0(w_bm(iwet))*IU;
            FLUX_bb = d0(w_bb(iwet))*IU;

            Gout.PFD_b  = DIV*FLUX_b;
            Gout.PFD_r  = DIV*FLUX_r;
            Gout.PFD_bm = DIV*FLUX_bm;
            Gout.PFD_bb = DIV*FLUX_bb;
            
        elseif( strcmp(Ptype,'PIC') )
            w_d = -r*MSK;
            FLUX_d = d0(w_d(iwet))*IU;
            Gout.PFD_d = DIV*FLUX_d;
            
        elseif ( strcmp(Ptype,'bSi') )
            w_d  = MSK;
            w_at = MSK;
            w_bt = MSK;
            
            w_at_tmp = -d*exp(-bt./T);
            w_bt_tmp = (d*r)./T;
            w_d(:,:,2:end-1)  = -r(:,:,2:end);
            w_at(:,:,2:end-1) = w_at_tmp(:,:,2:end);
            w_bt(:,:,2:end-1) = w_bt_tmp(:,:,2:end);
            
            FLUX_d  = d0(w_d(iwet))*IU;
            FLUX_at = d0(w_at(iwet))*IU;
            FLUX_bt = d0(w_bt(iwet))*IU;
            
            Gout.PFD_d  = DIV*FLUX_d;
            Gout.PFD_at = DIV*FLUX_at;
            Gout.PFD_bt = DIV*FLUX_bt;
        end
    end

    if (nargout > 2) % compute the hessian w.r.t. the parameters
        if ( strcmp(Ptype,'POP') | strcmp(Ptype,'POC') )
            a_b_b = 2*r./(b.^3);           % derivative of a_b w.r.t. b
            a_r_r = 0;                     % derivative of a_r w.r.t. r
            a_b_r = -1./b.^2;              % derivative of a_b w.r.t r
            a_bm_bm = (2*r.*T.^2)./(b.^3); % derivative of a_bm w.r.t. bm
            a_bb_bb = 2*r./b.^3;     % derivative of a_bb w.r.t. bb
            a_bm_bb = 2*r*T./b.^3;   % derivative of a_bm w.r.t. bb
            w_b_b = -a_b_b.*M;       % derivative of w_b w.r.t. b
            w_r_r = -a_r_r.*M;       % derivative of w_r w.r.t. r
            w_b_r = -a_b_r.*M;       % derivative of w_b w.r.t. r
            w_bm_bm = -a_bm_bm.*M;   % derivative of w_bm w.r.t. bm
            w_bb_bb = -a_bb_bb.*M;   % derivative of w_bb w.r.t. bb
            w_bm_bb = -a_bm_bb.*M;   % derivative of w_bm w.r.t. bb

            FLUX_b_b = d0(w_b_b(iwet))*IU;
            FLUX_bm_bm = d0(w_bm_bm(iwet))*IU;
            FLUX_bb_bb = d0(w_bb_bb(iwet))*IU;
            FLUX_b_r = d0(w_b_r(iwet))*IU;
            FLUX_bm_bb = d0(w_bm_bb(iwet))*IU;
            FLUX_r_r = 0;

            Hout.PFD_b_b   = DIV*FLUX_b_b;
            Hout.PFD_r_r   = DIV*FLUX_r_r;
            Hout.PFD_b_r   = DIV*FLUX_b_r;
            Hout.PFD_bm_bm = DIV*FLUX_bm_bm;
            Hout.PFD_bb_bb = DIV*FLUX_bb_bb;
            Hout.PFD_bm_bb = DIV*FLUX_bm_bb;
            
        elseif ( strcmp(Ptype,'PIC') )
            w_d_d = 0;
            FLUX_d_d = w_d_d*IU;
            Hout.PFD_d_d = DIV*FLUX_d_d;
            
        elseif ( strcmp(Ptype,'bSi') )
            r_at_at =  0;
            r_bt_bt =  r./T.^2;
            r_at_bt = -exp(-bt./T)./T ;
            w_d_d   =  0;
            w_at_at = -d.*r_at_at ;
            w_bt_bt = MSK;
            w_at_bt = MSK;
            w_at_d  = MSK;
            w_bt_d  = MSK; 
            w_bt_bt(:,:,2:end-1) = -d.*r_bt_bt(:,:,2:end) ;
            w_at_bt(:,:,2:end-1) = -d.*r_at_bt(:,:,2:end) ;
            w_at_d_tmp = -exp(-bt./T) ;
            w_at_d(:,:,2:end-1)  = w_at_d_tmp(:,:,2:end) ;
            w_bt_d_tmp = r./T ;
            w_bt_d(:,:,2:end-1)  = w_bt_d_tmp(:,:,2:end) ; 

            FLUX_d_d = d0(w_d_d)*IU;
            FLUX_at_at = d0(w_at_at)*IU;
            FLUX_bt_bt = d0(w_bt_bt(iwet))*IU;
            FLUX_at_bt = d0(w_at_bt(iwet))*IU;
            FLUX_at_d  = d0(w_at_d(iwet))*IU;
            FLUX_bt_d  = d0(w_bt_d(iwet))*IU;

            Hout.PFD_d_d   = DIV*FLUX_d_d;
            Hout.PFD_at_at = DIV*FLUX_at_at;
            Hout.PFD_bt_bt = DIV*FLUX_bt_bt;
            Hout.PFD_at_bt = DIV*FLUX_at_bt;
            Hout.PFD_at_d  = DIV*FLUX_at_d;
            Hout.PFD_bt_d  = DIV*FLUX_bt_d; 
        end
    end

end

function [DIV,IU,M,MSK,iwet] = mkOperators(M3d,grd);
% add an exra layer of zeros at the bottom to ensure there is no
% flux in or out of the domain when using periodic shift operators
% for finite differences and averaging
    [ny,nx,nz] = size(M3d);
    M3D        = zeros(ny,nx,nz+1);
    M3D(:,:,1:end-1) = M3d;

    % add the zw coordinate at the top of the extra layer
    ZW3d = grd.ZW3d;
    ZW3d = ZW3d(:,:,[1:end,end]);
    ZW3d(:,:,end) = grd.ZW3d(:,:,end)+grd.dzt(end);

    % areas of the top of the grid box
    dAt = (grd.DXT3d.*grd.DYT3d).*M3d;
    % volume of the grid boxes
    dVt = (grd.DXT3d.*grd.DYT3d.*grd.DZT3d).*M3d;
    %
    n = nx*ny*(nz+1);
    I0 = speye(n);
    i0 = zeros(ny,nx,nz+1);
    i0(:) = 1:n;
    % periodic shifts OK because M3D has a layer of zeros on the bottom
    iu = i0(:,:,[nz+1,1:nz]); %use a periodic upward shift
    ib = i0(:,:,[2:nz+1,1]); % use a periodic downward shift
    IU = I0(iu,:);
    IB = I0(ib,:);
    % keep only wet boxes
    iwet = find(M3D(:));
    I0 = I0(iwet,:); I0 = I0(:,iwet);
    IU = IU(:,iwet); IU = IU(iwet,:);
    IB = IB(:,iwet); IB = IB(iwet,:);

    % (compute the divergence in the center of the grid boxes)
    DIV = d0(dVt(iwet))\(I0-IB)*d0(dAt(iwet));

    % make a mask that zero's out the flux in through the sea surface
    % and out through the bottom of the ocean
    MSK = M3D.*M3D(:,:,[nz+1,1:nz]);
    M   = MSK.*ZW3d;
end

