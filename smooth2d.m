function [R] = smooth2d(y,r,m,grd,M3d);
d0 = @(x) spdiags(x(:),0,length(x(:)),length(x(:)));
msk = squeeze(M3d(:,:,1));
II = 0*msk;
n = length(II(:));
II(:) = 1:n;
ie = II(:,[2:end,1]); ie = ie(:);
iw = II(:,[end,1:end-1]); iw = iw(:);
in = II([2:end,1],:); in = in(:);
is = II([end,1:end-1],:); is = is(:);
i0 = II; i0 = i0(:);
I0 = speye(n);
IE = I0(ie,:); IW = I0(iw,:); IN = I0(in,:); IS = I0(is,:);
uDXt = d0(1./grd.DXU)*(IE-I0);
ebc = d0(msk(i0).*msk(ie));
uDXt = ebc*uDXt;
DX2 = d0(1./grd.DXT)*(I0-IW)*uDXt;
vDYt = d0(1./grd.DYV)*(IN-I0);
nbc = d0(msk(i0).*msk(in));
vDYt = nbc*vDYt;
DY2 = d0(1./grd.DYT)*(I0-IS)*vDYt;
L = I0-(r^2)*(DX2+DY2);
iwet = find(msk(:));
L = L(iwet,:); L  = L(:,iwet);
FL = mfactor(L);
R = nan*y; R(iwet) = y(iwet);

for k = 1:m
    R(iwet) = mfactor(FL,R(iwet));
end
