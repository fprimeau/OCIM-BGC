function [G,Gx] = uptake(parm, par)
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
c2p = parm.c2p;
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
Gpx(:,par.pindx.lbeta) = beta*alpha.*c2p.*dLdbeta;   % dGdlog_beta

% Gradient
% grad DIP
DIPx = zeros(nwet,nx);
np = 0;
if (par.biogeochem.opt_sigma == on)
    np = np + 1;
    DIPx(:,par.pindx.lsigma) = parm.Px(1:nwet,par.pindx.lsigma);
end

if (par.biogeochem.opt_kappa_dp == on)
    np = np + 1;
    DIPx(:,par.pindx.lkappa_dp) = parm.Px(1:nwet,par.pindx.lkappa_dp);
end

if (par.biogeochem.opt_slopep == on)
    np = np + 1;
    DIPx(:,par.pindx.lslopep) = parm.Px(1:nwet,par.pindx.lslopep);
end

if (par.biogeochem.opt_interpp == on)
    np = np + 1;
    DIPx(:,par.pindx.linterpp) = parm.Px(1:nwet,par.pindx.linterpp);

end

if (par.biogeochem.opt_alpha == on)
    np = np + 1;
    DIPx(:,par.pindx.lalpha) = parm.Px(1:nwet,par.pindx.lalpha);
end

if (par.biogeochem.opt_beta == on)
    np = np + 1;
    DIPx(:,par.pindx.lbeta) = parm.Px(1:nwet,par.pindx.lbeta);   
end

Gx = zeros(nwet,nx);
Gx = d0(Gp)*DIPx+d0(DIP)*Gpx;

% Hessian
% map P-model parameter order to N-model parameter order

% n = (nx+1)*nx/2;
% pos = tril(ones(nx),0);
% pos(pos==1) = 1:n;
% pos = pos';

% grad grad DIP
% DIPxx = zeros(n_wet,n);
% nn = 1;
% for ji = 1:np
    % for kk = ji:np
        % map(nn) = pos(ji,kk);
        % nn = nn + 1;
    % end
% end

% map = [pos(1,1) pos(1,2) pos(1,3) pos(1,4) pos(2,2) pos(2,3) pos(2,4) pos(3,3) pos(3,4) pos(4,4)];

% DIPxx(:,map) = parm.Pxx(1:n_wet,:);

% Gpxx = zeros(n_wet,n);
% Gpxx(:,pos(3,3)) = alpha*c2p.*L;                          % alpha alpha
% Gpxx(:,pos(3,4)) = beta*alpha*c2p.*dLdbeta;               % alpha beta
% Gpxx(:,pos(3,9)) = c2p_l*alpha*dc2pdc2p_l.*L;            % alpha c2p_l 
% Gpxx(:,pos(3,10)) = S*alpha*dc2pdS.*L;                    % alpha S 
% Gpxx(:,pos(4,4)) = beta*alpha*c2p.*dLdbeta+...          
    % beta^2*alpha*c2p.*d2Ldbetadbeta;                      % beta beta 
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





