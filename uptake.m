function [G,Gx] = uptake(parm)
% unpack the parameters to be optimized
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
Gpx(:,4) = alpha.*c2p.*L;              % dGdlog_alpha
Gpx(:,5) = beta*alpha.*c2p.*dLdbeta;   % dGdlog_beta

% Gradient
% grad DIP
DIPx = zeros(nwet,nx);
DIPx(:,[1,3:5]) = parm.Px(1:nwet,:);   

Gx = zeros(nwet,nx);
Gx = d0(Gp)*DIPx+d0(DIP)*Gpx;

% Hessian
% map P-model parameter order to N-model parameter order

% nx = parm.nx;
% n = (nx+1)*nx/2;
% pos = tril(ones(nx),0);
% pos(pos==1) = 1:n;
% pos = pos';

% grad grad DIP
% DIPxx = zeros(n_wet,n);
% map = [pos(1,1) pos(1,2) pos(1,3) pos(1,4) pos(2,2) pos(2,3) pos(2,4) ...
       % pos(3,3) pos(3,4) pos(4,4)];
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





