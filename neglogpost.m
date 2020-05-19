function [f, fx] = neglogpost(x, parm)
nx = length(x); % number of parameters
%
%++++++++++ check if the optimization routine suggests strange
% parameter values
A = exist('x0');
if A == 0
    x0 = parm.x0;
end
% xold = [x0(1 : 2); exp(x0(3 : end))];
% xnew = [x(1 : 2); exp(x(3 : end))];
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
xp = x([1 , 3 : 5]);
ip = [1, 6, 7, 8]; % index corresponding to all parameters;
[P, Px, parm] = eqPcycle(parm, ip, xp);

parm.P  = P ;
parm.Px = Px;

DIP = M3d+nan;  DIP(iwet) = P(1+0*nwet:1*nwet) ;
POP = M3d+nan;  POP(iwet) = P(1+1*nwet:2*nwet) ;
DOP = M3d+nan;  DOP(iwet) = P(1+2*nwet:3*nwet) ;
parm.DIP = DIP(iwet);

%%%%%%%%%%   End Solve P    %%%%%%%%%%%

%%%%%%%%%%%  Solve C   %%%%%%%%%%%%%%%%
xc = x;
% index corresponding to all parameters;
ic = [1,3,6,7,8,9,10,11];
[C,Cx] = eqCcycle(xc,parm,ic);
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

px  = zeros(nwet, nx);
px(:, [1, 3 : 5]) = Px(1 : nwet, :);

cx = Cx(1 : nwet, :);

if (nargout > 1)
    fx = ep.'*W*px + r*(ec.'*W*cx);
end
