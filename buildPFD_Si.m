function buildPFD_Si = buildPFD_Si(par,parm,Ptype);
% this code is used to build Particle Flux Diverngence (PFD)
%
%                      ________________________                                        
%                     /         A             /|  POP sinking          
%                 top/_______________________/ |       |       
%                    |  |                    | |       | -w   
%                    | /        V            | /       |       
%                 bot|/______________________|/        V
%                                                  
% PFD = (A*w(top)*POP(top)-A*w(bot)*POP(bot))/V;
% add an exra layer of zeros at the bottom to ensure there is no
% flux in or out of the domain when using periodic shift operators
% for finite differences and averaging
M3d = parm.M3d;
grd = parm.grd;
[ny,nx,nz] = size(M3d);
M3D = zeros(ny,nx,nz+1);
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
% (averages POP onto the top of the grid boxes)
AVG = d0((I0+IU)*M3D(iwet))\(I0+IU);
% (compute the divergence in the center of the grid boxes)
DIV = d0(dVt(iwet))\(I0-IB)*d0(dAt(iwet));
% (compute the flux at the top of the grid cells)
% mimics a Martin curve flux attenuation profile
%(see Kriest and Oschelies 2008 in Biogeosciences)

if (strcmp(Ptype,'POP')|strcmp(Ptype,'POC'))
    Tobs       = parm.Tobs;
    aveT       = parm.aveT;
    slope  = varargin{1};
    interp = varargin{2};
    T = aveT;
    % use slope = 0 for constant b value
    b = slope*T+interp;
elseif (strcmp(Ptype,'bSi'))
    at  = parm.at;
    bt  = parm.bt;
    T  = parm.sst(iwet) + 273.15;
    kappa_si = at*exp(-bt./T); 
elseif (strcmp(Ptype,'PIC'))
end
% ++++++++++++++++++++++++++++++++++++++++++++++++++

% ++++++++++++++++++++++++++++++++++++++++++++++++++

r = kappa_si;
b = parm.bsi;
a = r./b;
% particle sinking velocity at the top of the grid cells.
MSK = M3D.*M3D(:,:,[nz+1,1:nz]);
M = MSK.*ZW3d;

w = -a.*M(iwet);
%
dadb = -r./(b.^2);
dadr = 1./b;
drdat = exp(-bt./T);
drdbt = -(at*exp(-bt./T))./T;
dadat =  dadr.*drdat;
dadbt =  dadr.*drdbt;
%
dwdb = -dadb.*M(iwet);
dwdat = -dadat.*M(iwet);
dwdbt = -dadbt.*M(iwet);
%
d2adb2 = 2*r./(b^3);
d2adat2 = 0;
d2adbt2 = dadr.*(at*exp(-bt./T))./(T.^2);
d2adatdbt = -dadr.*exp(-bt./T)./T;
d2adatdb  = -(1./b.^2).*drdat;
d2adbtdb  = -(1./b.^2).*drdbt;

d2wdb2 = -d2adb2.*M(iwet);
d2wdat2 = 0;
d2wdbt2 = -d2adbt2.*M(iwet);
d2wdatdbt = -d2adatdbt.*M(iwet);
d2wdatdb  = -d2adatdb.*M(iwet);
d2wdbtdb  = -d2adbtdb.*M(iwet);
%FLUX = d0(w(iwet))*AVG;
FLUX = d0(w)*IU;
dFLUXdb = d0(dwdb)*IU;
dFLUXdat = d0(dwdat)*IU;
dFLUXdbt = d0(dwdbt)*IU;
%
d2FLUXdb2 = d0(d2wdb2)*IU;
d2FLUXdat2 = d0(d2wdat2)*IU;
d2FLUXdbt2 = d0(d2wdbt2)*IU;
d2FLUXdatdbt = d0(d2wdatdbt)*IU;
d2FLUXdatdb  = d0(d2wdatdb)*IU;
d2FLUXdbtdb  = d0(d2wdbtdb)*IU;
% particle flux divergence operator
vout.PFdiv = DIV*FLUX;
vout.dPFDdb = DIV*dFLUXdb;
vout.dPFDdat = DIV*dFLUXdat;
vout.dPFDdbt = DIV*dFLUXdbt;
%
vout.d2PFDdb2 = DIV*d2FLUXdb2;
vout.d2PFDdat2 = DIV*d2FLUXdat2;
vout.d2PFDdbt2 = DIV*d2FLUXdbt2;
vout.d2PFDdatdbt = DIV*d2FLUXdatdbt;
vout.d2PFDdatdb = DIV*d2FLUXdatdb;
vout.d2PFDdbtdb = DIV*d2FLUXdbtdb;
%% +++++++++++++++++++++++++++++++++++++++++++