function [PFD,dPFDdd,d2PFDdd2] = buildPFD_CaCO3(p,grd,M3d,iwet);
% Particle Flux Divergence operator
% compute a sinking velocity that produces an exponentially
% decaying PIC flux
iwet = find(M3d(:));
n=length(iwet);
[ny,nx,nz] = size(M3d);
d0 = @(v) spdiags([v(:)],[0],length(v(:)),length(v(:)));
%
I0 = speye(n);
i0 = zeros(ny,nx,nz);
i0(iwet) = 1:length(iwet); % each ocean point on the grd gets unique index  
iu = i0(:,:,[1,1:nz-1]); % index of neighbor above
msk = M3d;
msk(:,:,1) = 0; % impose no flux through the top
IU = d0(msk(iwet))*I0(iu(iwet),:);

% 
dVt = grd.DXT3d.*grd.DYT3d.*grd.DZT3d;
D = grd.ZT3d.*M3d;
V = dVt.*M3d;
A = (grd.DXT3d.*grd.DYT3d).*M3d;
d = p.d; % vertical e-folding length scale for PIC flux attenuation
r = 1./p.taup;
w = -r*d;
dwdd = -r;
d2wdd2 = 0;
%
msk = M3d(:,:,[2:end,end]);
msk(:,:,end) = 0;
%
fb = (msk.*A).*w;
dfbdd = (msk.*A).*dwdd;
d2fbdd2 = (msk.*A).*d2wdd2;
PFD = d0(V(iwet))\((IU-I0)*d0(fb(iwet)));
dPFDdd = d0(V(iwet))\((IU-I0)*d0(dfbdd(iwet)));
d2PFDdd2 = d0(V(iwet))\((IU-I0)*d0(d2fbdd2(iwet)));