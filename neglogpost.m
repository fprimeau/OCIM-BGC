function [f, fx] = neglogpost(x, parm, par)
on = true; off = false;
nx = length(x); % number of parameters
%
%++++++++++ check if the optimization routine suggests strange
% parameter values
A = exist('x0');
if A == 0
    x0 = parm.p0;
end

xold = [exp(x0(1 : end))];
xnew = [exp(x(1 : end))];
%
ibad = find(xnew > 5*xold | xnew < 1/5*xold);
x(ibad) = x0(ibad);
x0 = x;
%+++++++++restore the parameter values back to their original ones.
xnew = [exp(x(1 : end))];
fprintf('current parameter is:  \n');
for ii = 1 : length(x)
    fprintf('%3.3e;  ',xnew(ii));
end
fprintf('\n')

dVt  = parm.dVt  ;
M3d  = parm.M3d  ;
iwet = parm.iwet ;
nwet = parm.nwet ;

%%%%%%%%%%%%   Slove P    %%%%%%%%%%%%%
[P, Px, parm] = eqPcycle(par, parm, x);

parm.P  = P ;
parm.Px = Px;

DIP = M3d+nan;  DIP(iwet) = P(1+0*nwet:1*nwet) ;
POP = M3d+nan;  POP(iwet) = P(1+1*nwet:2*nwet) ;
DOP = M3d+nan;  DOP(iwet) = P(1+2*nwet:3*nwet) ;
parm.DIP = DIP(iwet);

%%%%%%%%%%   End Solve P    %%%%%%%%%%%

%%%%%%%%%%%  Solve C   %%%%%%%%%%%%%%%%
[C,Cx] = eqCcycle(par, parm, x);
DIC = M3d+nan; DIC(iwet) = C(0*nwet+1:1*nwet) ;
POC = M3d+nan; POC(iwet) = C(1*nwet+1:2*nwet) ;
ODC = M3d+nan; DOC(iwet) = C(2*nwet+1:3*nwet) ;
CaC = M3d+nan; CaC(iwet) = C(3*nwet+1:4*nwet) ;
%%%%%%%%%%%%% End solve C %%%%%%%%%%%%%%%

W   = d0(dVt(iwet)/sum(dVt(iwet)));

DIC = DIC + parm.human_co2;
ep  = DIP(iwet) - parm.po4obs(iwet);
ec  = DIC(iwet) - parm.DICobs(iwet);

r   = 0.00093932848354556/99.382557451999588; % from Teng et al.(2014);
f   = 0.5*(ep.'*W*ep) + 0.5*r*(ec.'*W*ec);

fprintf('current objective function value is %3.3e \n',f);

if (nargout > 1)
    if (par.biogeochem.opt_sigma == on)
        px(:, par.pindx.lsigma) = Px(1:nwet,par.pindx.lsigma);
        cx(:, par.pindx.lsigma) = Cx(1:nwet,par.pindx.lsigma);
    end

    if (par.biogeochem.opt_slopep == on)
        px(:, par.pindx.lslopep) = Px(1:nwet,par.pindx.lslopep);
        cx(:, par.pindx.lslopep) = Cx(1:nwet,par.pindx.lslopep);
    end
    
    if (par.biogeochem.opt_interpp == on)
        px(:, par.pindx.linterpp) = Px(1:nwet,par.pindx.linterpp);
        cx(:, par.pindx.linterpp) = Cx(1:nwet,par.pindx.linterpp);
    end
    
    if (par.biogeochem.opt_alpha == on)
        px(:, par.pindx.lalpha) = Px(1:nwet,par.pindx.lalpha);
        cx(:, par.pindx.lalpha) = Cx(1:nwet,par.pindx.lalpha);
    end
    
    if (par.biogeochem.opt_beta == on)
        px(:, par.pindx.lbeta) = Px(1:nwet,par.pindx.lbeta);
        cx(:, par.pindx.lbeta) = Cx(1:nwet,par.pindx.lbeta);
    end
    
    if (par.biogeochem.opt_kappa_dp == on)
        px(:, par.pindx.lkappa_dp) = Px(1:nwet,par.pindx.lkappa_dp);
        cx(:, par.pindx.lkappa_dp) = Cx(1:nwet,par.pindx.lkappa_dp);
    end
    
    if (par.biogeochem.opt_slopec == on)
        px(:, par.pindx.lslopec) = zeros(nwet,1);
        cx(:, par.pindx.lslopec) = Cx(1:nwet,par.pindx.lslopec);
    end
    
    if (par.biogeochem.opt_interpc == on)
        px(:, par.pindx.linterpc) = zeros(nwet,1);
        cx(:, par.pindx.linterpc) = Cx(1:nwet,par.pindx.linterpc);
    end
    
    if (par.biogeochem.opt_d == on)
        px(:, par.pindx.ld) = zeros(nwet,1);
        cx(:, par.pindx.ld) = Cx(1:nwet,par.pindx.ld);
    end
    
    if (par.biogeochem.opt_kappa_dc == on)
        px(:, par.pindx.lkappa_dc) = zeros(nwet,1);
        cx(:, par.pindx.lkappa_dc) = Cx(1:nwet,par.pindx.lkappa_dc);
    end
    
    if (par.biogeochem.opt_RR == on)
        px(:, par.pindx.lRR) = zeros(nwet,1);
        cx(:, par.pindx.lRR) = Cx(1:nwet,par.pindx.lRR);
    end
    
    fx = ep.'*W*px + r*(ec.'*W*cx);
end
