function SM = smoothit(grd,M3d,data,n,r)
% this function is used to smooth 2d field
% grd is model grid;
% M3d is ocean-land masks
% data is the 2d field to be smoothed
% n is a integer number
% r is a smooth length scale. The larger n and r are, the more smooth.
tmp = M3d(:,:,1);
iwet = find(tmp(:)==1);

y = data(iwet);
y = inpaint_nans(y);
data(iwet) = y;

SM = smooth2d(data,r,n,grd,M3d); % increase r and/or n for soother field



