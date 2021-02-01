%% simplifying cell params hessian like for P 
% P model only
for jj = 1:npx
            for jk = jj:npx
                kk = kk + 1;
                if (par.opt_sigma == on)
                    % sigma sigma
                    if (jj == jk & jj == pindx.lsigma)
                        tmp = sigma*[RR*G*C2P + 2*d0(RR*Gx(:,jj))*C2P; ...
                                     -G*C2P - 2*d0(Gx(:,jj))*C2P;...
                                      G*C2P + 2*d0(Gx(:,jj))*C2P;...
                                     -RR*G*C2P - 2*d0(RR*Gx(:,jj))*C2P; ...
                                     2*RR*G*C2P + 2*2*d0(RR*Gx(:,jj))*C2P] + ...
                              [-d0((I+(1-sigma)*RR)*Gxx(:,kk))*C2P; ...
                               (1-sigma)*d0(Gxx(:,kk))*C2P; ...
                               sigma*d0(Gxx(:,kk))*C2P; ...
                               d0((1-sigma)*RR*Gxx(:,kk))*C2P; ...
                               d0(-2*(1-sigma)*RR*Gxx(:,kk))*C2P + ...
                               d0(N2C*Gxx(:,kk))*C2P];

                        Cxx(:,kk) = mfactor(FD, tmp);
                        %pairs not assciated with sigma;
                    elseif (jj ~= pindx.lsigma & jk ~= pindx.lsigma)
                        tmp = [-d0((I+(1-sigma)*RR)*Gxx(:,kk))*C2P; ...
                               (1-sigma)*d0(Gxx(:,kk))*C2P; ...
                               sigma*d0(Gxx(:,kk))*C2P; ...
                               d0((1-sigma)*RR*Gxx(:,kk))*C2P; ...
                               d0(-2*(1-sigma)*RR*Gxx(:,kk))*C2P + ...
                               d0(N2C*Gxx(:,kk))*C2P];

                        Cxx(:,kk) = mfactor(FD, tmp);

                    else
                        tmp = [RR*sigma*d0(Gx(:,jk))*C2P; ...
                               -sigma*d0(Gx(:,jk))*C2P; ...
                               sigma*d0(Gx(:,jk))*C2P; ...
                               -RR*sigma*d0(Gx(:,jk))*C2P; ...
                               2*RR*sigma*d0(Gx(:,jk))*C2P] + ...
                              [-d0((I+(1-sigma)*RR)*Gxx(:,kk))*C2P; ...
                               (1-sigma)*d0(Gxx(:,kk))*C2P; ...
                               sigma*d0(Gxx(:,kk))*C2P; ...
                               d0((1-sigma)*RR*Gxx(:,kk))*C2P;...
                               d0(-2*(1-sigma)*RR*Gxx(:,kk))*C2P + ...
                               d0(N2C*Gxx(:,kk))*C2P];

                        Cxx(:,kk) = mfactor(FD, tmp);
                    end
                else
                    tmp = [-d0((I+(1-sigma)*RR)*Gxx(:,kk))*C2P; ...
                           (1-sigma)*d0(Gxx(:,kk))*C2P; ...
                           sigma*d0(Gxx(:,kk))*C2P; ...
                           d0((1-sigma)*RR*Gxx(:,kk))*C2P; ...
                           d0(-2*(1-sigma)*RR*Gxx(:,kk))*C2P + ...
                           d0(N2C*Gxx(:,kk))*C2P];

                    Cxx(:,kk) = mfactor(FD, tmp);
                    % sigma foo
                end
            end
end

% Q10Photo Q10Photo
        if (par.opt_Q10Photo)
			C2P_Q10_Q10 = par.CellOut.d2C2P_dQ10Photo2(iwet);
            kk = kk + 1;
            tmp = par.BIO.Q10Photo*[-(I+(1-sigma)*RR)*(G*(C2P_Q10+par.BIO.Q10Photo*C2P_Q10_Q10)); ...
                      (1-sigma)*G*(C2P_Q10+par.BIO.Q10Photo*C2P_Q10_Q10); ...
                      sigma*G*(C2P_Q10+par.BIO.Q10Photo*C2P_Q10_Q10); ...
                      (1-sigma)*RR*(G*(C2P_Q10+par.BIO.Q10Photo*C2P_Q10_Q10)); ...
                      -2*(1-sigma)*RR*(G*(C2P_Q10+par.BIO.Q10Photo*C2P_Q10_Q10)) + ...
                      N2C*G*(C2P_Q10+par.BIO.Q10Photo*C2P_Q10_Q10)];

            Cxx(:,kk) = mfactor(FD, tmp);
        end
% Q10Photo Q10Photo
        if (par.opt_Q10Photo)
			C2P_Q10_Q10 = par.CellOut.d2C2P_dQ10Photo2(iwet);
            dC2Ptmp = C2P_Q10+par.BIO.Q10Photo*C2P_Q10_Q10;
            kk = kk + 1;
            tmp = par.BIO.Q10Photo*[-(I+(1-sigma)*RR)*(G*(dC2Ptmp)); ...
                      (1-sigma)*G*(dC2Ptmp); ...
                      sigma*G*(dC2Ptmp); ...
                      (1-sigma)*RR*(G*(dC2Ptmp)); ...
                      -2*(1-sigma)*RR*(G*(dC2Ptmp)) + ...
                      N2C*G*(dC2Ptmp)];

            Cxx(:,kk) = mfactor(FD, tmp);
        end
        
        % Q10Photo fStorage
        if (par.opt_Q10Photo & par.opt_fStorage)
			C2P_Q10_fStor = par.CellOut.d2C2P_dQ10Photo_dfStorage(iwet);
            dC2Ptmp=C2P_Q10_fStor;
            kk = kk + 1;
            tmp = par.BIO.Q10Photo*par.BIO.fStorage*[-(I+(1-sigma)*RR)*(G*dC2Ptmp); ...
                         (1-sigma)*G*dC2Ptmp; ...
                         sigma*G*dC2Ptmp; ...
                         (1-sigma)*RR*(G*dC2Ptmp); ...
                         -2*(1-sigma)*RR*(G*dC2Ptmp) + ...
                         N2C*G*dC2Ptmp];

            Cxx(:,kk) = mfactor(FD, tmp);
        end
        % % Q10Photo fStorage
        % if (par.opt_Q10Photo & par.opt_fStorage)
		% 	C2P_Q10_fStor = par.CellOut.d2C2P_dQ10Photo_dfStorage(iwet);
        %     kk = kk + 1;
        %     tmp = par.BIO.Q10Photo*par.BIO.fStorage*[-(I+(1-sigma)*RR)*(G*C2P_Q10_fStor); ...
        %                  (1-sigma)*G*C2P_Q10_fStor; ...
        %                  sigma*G*C2P_Q10_fStor; ...
        %                  (1-sigma)*RR*(G*C2P_Q10_fStor); ...
        %                  -2*(1-sigma)*RR*(G*C2P_Q10_fStor) + ...
        %                  N2C*G*C2P_Q10_fStor];
		% 
        %     Cxx(:,kk) = mfactor(FD, tmp);	
        % end

        % PLip_PCutoff PStor_rCutoff
        if (par.opt_PLip_PCutoff & par.opt_PStor_rCutoff)
			C2P_PCutoff_rCutoff = par.CellOut.d2C2P_drCutoff_dPCutoff(iwet);
            kk = kk + 1;
            tmp = par.BIO.PLip_PCutoff*par.BIO.PStor_rCutoff*[-(I+(1-sigma)*RR)*(G*C2P_PCutoff_rCutoff); ...
                         (1-sigma)*G*C2P_PCutoff_rCutoff; ...
                         sigma*G*C2P_PCutoff_rCutoff; ...
                         (1-sigma)*RR*(G*C2P_PCutoff_rCutoff); ...
                         -2*(1-sigma)*RR*(G*C2P_PCutoff_rCutoff) + ...
                         N2C*G*C2P_PCutoff_rCutoff];

            Cxx(:,kk) = mfactor(FD, tmp);
        end
        
%% why is fsolve  not working? 
  i=152
    c3 = [6.06597219280557e-05]
    c2 = [3.89440203334762e-06]
    c1= [6.03668036945044e-07]
    c0 = [-6.42627788907703e-08]
    rhsFuncColim = @(E) (c3*E.^3 + c2*E.^2 + c1*E +c0)
    
    E = [0:0.00001:0.1];
    y = rhsFuncColim(E);
    
    figure;
    plot(E,y)
    grid on
    xlabel('E'); ylabel('rhsFuncColim(E)')
    title('i=152')
    % zero at e= 0.060535
    
%%    
c3 = [6.06597219280557e-05];
    c2 = [3.89440203334762e-06];
    c1= [6.03668036945044e-07];
    c0 = [-6.42627788907703e-08];
    rhsFuncColim = @(E) (c3*E.^3 + c2*E.^2 + c1*E +c0)
    
    E0 = 0.22;
    fi = rhsFuncColim(E0);
    
    f = @(E) rhsFuncColim(E)/fi;
    
options = optimoptions('fsolve','Display','iter','FunctionTolerance',1e-9,'OptimalityTolerance',1e-9);
%[x fval] = fsolve(rhsFuncColim,0.01,options)
[x fval] = fsolve(f,E0,options)

%f(x) - f(x+1) <tol

%rhsFuncColim(0.2)-rhsFuncColim(0.19)

%%
bads = M3dsurf*0;
bads(iprod(ibad)) = 1;
%all wet points = iprod
%land pints = idry
bads(iprod)=bads(iprod)+1;
figure;
pcolor(1:2:360,-89:2:89, bads(:,:,1)); colorbar
title('fProtA>1 =bad=2, wet=1, dry = 0')


%% testing and correcting light field
spd  = 24*60^2 ;
PARobs=load('/Users/megansullivan/Documents/UC Irvine/DATASETS/weilei_gp_DATA/annual_PAR_90x180.mat')
PARobs_PPFD = PARobs.par*10^6/spd; % PAR at surface
clear PARobs


iwet1 = find(M3dsurf(:,:,1));
%ibad = find(isnan(PARobs_PPFD(iwet1)));

PARsurf        = 0*M3dsurf(:,:,1);
PARsurf(iwet1)=PARobs_PPFD(iwet1);

[ibady,ibadx] =find(isnan(PARsurf));

for ii=1:length(ibadx)
    ilat = ibady(ii);
    ilon = ibadx(ii);
    if ilon>1 & ilon <180
        nearby=PARsurf(ilat-1:ilat+1,ilon-1:ilon+1); % would be better to weight values on same latitude more
        
        if length(nearby(nearby>0)) > 0
            PARsurf(ilat,ilon) = nanmean(nearby(nearby>0));
        else
            nearby=PARsurf(ilat-2:ilat+2,ilon-2:ilon+2);
            if length(nearby(nearby>0)) > 0
                PARsurf(ilat,ilon) = nanmean(nearby(nearby>0));
            else
                keyboard
            end
        end
        
    else % if ilon==180
        nearby=PARsurf(ilat-1:ilat+1,ilon-2:ilon);
        if length(nearby(nearby>0)) > 0
            PARsurf(ilat,ilon) = nanmean(nearby(nearby>0));
        else
            keyboard
        end

    end
end


%     if PARgrd(ibady(ii),ibadx(ii)+1,1)>0
%         PARgrd(ibady(ii),ibadx(ii),1) = PARgrd(ibady(ii),ibadx(ii)+1,1);
%     elseif PARgrd(ibady(ii),ibadx(ii)-1,1)>0
%         PARgrd(ibady(ii),ibadx(ii),1) = PARgrd(ibady(ii),ibadx(ii)-1,1);
%     elseif PARgrd(ibady(ii)+1,ibadx(ii),1)>0
%         PARgrd(ibady(ii),ibadx(ii),1) = PARgrd(ibady(ii)+1,ibadx(ii),1);
%     elseif PARgrd(ibady(ii)-1,ibadx(ii),1)>0
%         PARgrd(ibady(ii),ibadx(ii),1) = PARgrd(ibady(ii)-1,ibadx(ii),1);
%     else
%         igood =find(PARgrd(ibady(ii),:,1));
%         PARgrd(ibady(ii),ibadx(ii),1) = nanmean(PARgrd(igood));
%     end
    

%par.kI = 0.04;   % Light attenuation coefficient in seawater [m^-1]
%	PAR        = 0*M3d;
%	for ii=1:par.nzo % only needed in euphotic zone for cell growth
%		PAR(:,:,ii) = PARobs_PPFD.*exp(-par.kI*grd.zt(ii)); %PAR at mid depth of grid box [ii]
%    end
%    par.PARobs = PAR;
%		clear PAR


%%
par.PARobs = M3dsurf*0;
par.PARobs(:,:,1)=PARsurf;
par.PARobs(:,:,2)=0.66*PARsurf;

Irr0 = par.PARobs(iprod);

%%
% check is same size as surface nitrate/phosphate. (nans/coastlines) -
% where nans don't match up
% is surface data same as model output (goodindx
%Irr0, N0, P0,
%iprod = find(M3d(:,:,1:2)); %production in top two layers
iprod= find(M3dsurf(:));
idry= find(~M3dsurf(:));

ibadIrr = find(isnan(Irr0));

ibadN = find(isnan(N0));

ibadP = find(isnan(P0));

ibadT = find(isnan(T0));

% problem is indeed the light field (288 NaNs in wet points)

isnan(Irr0)



%%
outPath='/Users/megansullivan/Documents/UC Irvine/DATASETS/weilei_gp_DATA/'


bads = M3dsurf*0;
bads(iprod(ibadIrr)) = 1;
%all wet points = iprod
%land pints = idry
bads(iprod)=bads(iprod)+1;
figure;
pcolor(1:2:360,-89:2:89, bads(:,:,1)); colorbar
title('annual_ PAR_ 90x180.mat PAR field: nan = 2, wet=1, dry = 0')


figTitle = ['annual_PAR_field_nans'];
%print(gcf,[outPath 'FIG_' figTitle '.png'],'-dpng')

%% 
PARgrd = M3dsurf*0;  % NaN(size(M3dsurf));
PARgrd(iprod)=Irr0;

[iy1,ix1]=find(isnan(PARgrd(:,:,1)));

PARgrd(iy1(1),ix1(1),1) = PARgrd(iy1(1),ix1(1)+1,1);

for ii=1:length(ix1)
    if PARgrd(iy1(ii),ix1(ii)+1,1)>0
        PARgrd(iy1(ii),ix1(ii),1) = PARgrd(iy1(ii),ix1(ii)+1,1);
    elseif PARgrd(iy1(ii),ix1(ii)-1,1)>0
        PARgrd(iy1(ii),ix1(ii),1) = PARgrd(iy1(ii),ix1(ii)-1,1);
    elseif PARgrd(iy1(ii)+1,ix1(ii),1)>0
        PARgrd(iy1(ii),ix1(ii),1) = PARgrd(iy1(ii)+1,ix1(ii),1);
    elseif PARgrd(iy1(ii)-1,ix1(ii),1)>0
        PARgrd(iy1(ii),ix1(ii),1) = PARgrd(iy1(ii)-1,ix1(ii),1);
    else
        igood =find(PARgrd(iy1(ii),:,1));
        PARgrd(iy1(ii),ix1(ii),1) = nanmean(PARgrd(igood));
    end
    
end

find(isnan(PARgrd(:,:,1)))
% now interpolate light field to bads in iprod 
% first look at plot of light. make sure straight extrapolation on lines of
% latitude will be close enough

%[ix,iy]=find(bads(:,:,1)==2);

% need Irr0 in matrix form (before taking iprod indexes)
%Irr0(ibadIrr)


%% load inputs
load('/Users/megansullivan/Documents/UC Irvine/DATASETS/TraitModel_TestData/inputs_surf_CellCNP.mat')
par.BIO = parBIO;

on = true; off = false;
par.opt_Q10Photo = on;
par.opt_fRibE 	 = off;
par.opt_fStorage = off;
par.opt_kST0 	 = off;

par.pindx = pindx;
%% plot no3

no3in = M3dsurf*NaN;
no3in(iprod) = N0;

figure; hold on
contourf(1:2:360,-89:2:89,no3in(:,:,1)); cb=colorbar;
colormap(flipud(autumn))

title('Nitrate input field','Fontsize',18);
xlabel('Longitude');
ylabel('Latitude');
ylabel(cb,'NO3 [g/m^3]');
grid off

figTitle = 'NO3_surface';
print(gcf,[outPath 'FIG_' figTitle '_' ver '.png'],'-dpng')


%% P0
po4in = M3dsurf*NaN;
po4in(iprod) = P0;

figure; hold on
contourf(1:2:360,-89:2:89,po4in(:,:,1)); cb=colorbar;
colormap(flipud(autumn))

title('surface P input field','Fontsize',18);
xlabel('Longitude');
ylabel('Latitude');
ylabel(cb,'PO4 [g/m^3]');
grid off

figTitle = 'PO4_surface';
print(gcf,[outPath 'FIG_' figTitle '_' ver '.png'],'-dpng')