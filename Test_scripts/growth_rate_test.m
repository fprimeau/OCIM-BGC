% plot model vs. satellite

%% Load satellite growth rate
fname = '/Users/megansullivan/Documents/UC Irvine/DATASETS/TFM20_Monthly.nc';
finfo = ncinfo(fname)

%'MU_CBPM' monthly
mu_cbpm_monthly = ncread(fname,'MU_CBPM');
mu_cbpm_tmp = nanmean(mu_cbpm_monthly,4);


lon_cbpm = ncread(fname,'Longitude');
lat_cbpm = ncread(fname,'Latitude');

mu_grd = interpTo90x180(lon_cbpm,lat_cbpm,mu_cbpm_tmp,'FT');
mu_cbpm = mu_grd.m /24; %convert units for per day to per hour

load('/Users/megansullivan/Documents/UC Irvine/GitHub/TraitModel_output/FIGS_02May21/Int_CNPP.mat')

load('/Users/megansullivan/Documents/UC Irvine/GitHub/TraitModel_output/FIGS_02May21/mu_model.mat')

%iprod = find(mu_model~=0);
iwet = iprod;

%% plot 
mu_model1 = mu_model(:,:,1);

iprod = find(mu_model(:,:,1));

figure;
plot(mu_cbpm(iprod),mu_model1(iprod),'ok')
xlabel('CbPM mu [1/hr]');
ylabel('Model mu [1/hr]');
title('model growth rate comparison to CbPM')
xlim([0 0.07])
ylim([0 0.07])
grid on

%%
    %mu_cbpm
    imu   = find(mu_cbpm(:) > 0 & DOP(:)>0) ;
    iMU_ATL = find(mu_cbpm(iwet)>0 & ATL(iwet)>0) ;
    iMU_PAC = find(mu_cbpm(iwet)>0 & PAC(iwet)>0) ;
    iMU_IND = find(mu_cbpm(iwet)>0 & IND(iwet)>0) ;
    iMU_ARC = find(mu_cbpm(iwet)>0 & ARC(iwet)>0) ;
    iMU_MED = find(mu_cbpm(iwet)>0 & MED(iwet)>0) ;
    fprintf('R^2 for DOP is %3.3f \n',rsquare(mu_cbpm(imu),DOP(imu)))

    nfig = nfig + 1;
    figure(nfig)
    plot(mu_cbpm(iwet(iMU_ATL)), DOP(iwet(iMU_ATL)),'ro')
    hold on
    plot(mu_cbpm(iwet(iMU_PAC)), DOP(iwet(iMU_PAC)),'ks')
    hold on
    plot(mu_cbpm(iwet(iMU_IND)), DOP(iwet(iMU_IND)),'b^')
    hold on
    plot(mu_cbpm(iwet(iMU_ARC)), DOP(iwet(iMU_ARC)),'g*')
    hold on
    plot(mu_cbpm(iwet(iMU_MED)), DOP(iwet(iMU_MED)),'c>')
    hold on
	legend('ATL','PAC','IND','ARC','Location','northwest')
	hold on
    plot([0 0.75],[0 0.75],'r-','linewidth',3)
    xlim([0 0.75])
    ylim([0 0.75])

	xlabel('Observed DOP (mmol/m^3)');
    ylabel('Model DOP (mmol/m^3)');
	figTitle = 'DOPcompare2obs';
	print(gcf,[figPath 'FIG_' figTitle '.png'],'-dpng')

%%
function [out] =interpTo90x180(lon,lat,o2,method)
% this function can interpolate irregular CMIP5 models to
% [out]=interpTo360x180(lon,lat,o2,method)
% regular 1 degree x 1 degree grid
% lon and lat are model coordinate if 2D, must equal the size
% of o2, o2 could be 3D or 2D, the last dimensions are the
% same as specified in lon and lat
% method can be 'nearest' or 'FT'
% Detailed explanation goes here
% out is struct containing X, Y and interperlated matrix (m)
% out.x out.y, out.m


%%% method 1
if strcmp(method,'FT')

    %    [X,Y] = meshgrid(0:2:360,-90:2:90);
         [X,Y] = meshgrid(1:2:359,-89:2:89);
    %    [X,Y] = meshgrid(0.5:1:359.5,-89.5:1:89.5);
    %    [X,Y] = meshgrid(0.5:2:359.5,-89.5:2:89.5);
    out.x=X;out.y=Y;

    % squeeze o2 to get real dimensions
    o2=squeeze(o2);

    if ndims(o2) == 3
        %    temp=squeeze(o2(1,ind,:,:));
        %elseif ndims(o2) == 3
        %    temp=squeeze(o2(ind,:,:));
    else
        temp=squeeze(o2(:,:));
    end

    if size(lat,2)==1 | size(lat,1) == 1

        [X1,Y1] = meshgrid(wrapTo360(lon),lat);

        FT =TriScatteredInterp(X1(:), Y1(:), temp(:));

        out.m=FT(X,Y);
    else

        FT = TriScatteredInterp(wrapTo360(lon(:)),lat(:),temp(:));

        %index=find(isfinite(FT)==1);
        %mean=nansum(FT(index)*model.lat{ii}(index))/nansum(model.lat{ii}(index));

        out.m=FT(X,Y);
    end

    %get the mean of FT by accounting for latitude-dependence
end

%%% method 2
if strcmp(method,'nearest')

    %     X=0.5:1:359.5; Y=-89.5:1:89.5;
    X = 1:2:359; Y = -89:2:89; 
    out.x=X;out.y=Y;

     if ndims(o2) == 4
         temp=squeeze(o2(1,ind,:,:));
    elseif ndims(o2) == 3
        temp=squeeze(o2(ind,:,:));
    else
        temp=squeeze(o2(:,:));
    end

    for i=1:length(X)
        for j=1:length(Y)

            dist=abs(X(i)-wrapTo360(lon)) + abs(Y(j)-lat);

            dist_min=min(min(dist));

            [ix, jx]=find(dist==dist_min);

            out.m(j,i)=temp(min(ix),min(jx));

            %whos ix jx x;
            %temp(i,j) =x;
       end
   end
end
end