function [kw,p] = kw(M3d,grd);
    fname = 'gasx_ocmip2.nc';
    lon = ncread(fname,'LON')-180;
    lon = [lon(1)-1;lon(:);lon(end)+1];
    lat = ncread(fname,'LAT');
    [LON,LAT] = meshgrid(lon,lat);
    msk = ncread(fname,'TMASK');
    msk = msk';
    msk = [msk(:,181:end),msk(:,1:180)];
    [nlat,nlon] = size(msk);
    
    xkw = zeros(nlat,nlon,12);
    fice = zeros(nlat,nlon,12);
    p = zeros(nlat,nlon,12);
    
    tmp1 = ncread(fname,'FICE');
    tmp2 = ncread(fname,'XKW');
    tmp3 = ncread(fname,'P');
    for i = 1:12
        fice(:,:,i) = tmp1(:,:,i)';
        xkw(:,:,i) = tmp2(:,:,i)';
        p(:,:,i) = tmp3(:,:,i)';
    end
    fice = mean(fice,3); % average over months
    fice = [fice(:,181:end),fice(:,1:180)];
    xkw = mean(xkw,3);
    xkw = [xkw(:,181:end),xkw(:,1:180)];
    p = mean(p,3);
    p = [p(:,181:end),p(:,1:180)];
    
    kw = (1-fice).*xkw;
    ilnd = find(msk(:)==0);
    kw(ilnd) = NaN;
    
    p = [p(:,end),p,p(:,1)];
    kw = [kw(:,end),kw,kw(:,1)];
    kw = inpaint_nans(kw);
    % interpolate onto the model grd
    p = interp2(LON,LAT,p,grd.XT3d(:,:,1),grd.YT3d(:,:,1));
    kw = interp2(LON,LAT,kw,grd.XT3d(:,:,1),grd.YT3d(:,:,1));
    msk = M3d(:,:,1);
    ilnd = find(msk(:)==0);
    kw(ilnd) = NaN;
    