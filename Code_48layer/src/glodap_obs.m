function data = Gop(varargin)
% function data = Gop(data,vn,t0,tf)
%     data: required, struct  
%       vn: required, variable name for extract
%       t0: optional, starting year for data extraction
%       tf: optional, ending year for data extraction
% extract variable from GLODAPv2.2021 data
% interpolate data to OCIM grid
% build an operator to map OCIM data to where
% GLODAP data available

  p = inputParser;
  % addRequired(p,'data',@isnumeric);
  addRequired(p,'data',@isstruct);
  addRequired(p,'vn',@ischar);
  % addOptional(p,'t0',1850,@isstruct)
  addOptional(p,'t0',1750,@isnumeric)
  % addOptional(p,'tf',2020,@isstruct)
  addOptional(p,'tf',2022,@isnumeric)
  parse(p,varargin{:});
  p.Results

  fr   = 0;
  vn = p.Results.vn; % vn   = 'G2tco2';
  vnv  = sprintf('%sf',vn);
  vnw  = vn(3:end);
  parms.t0 = p.Results.t0;
  parms.tf = p.Results.tf;
  %load ~/Research/CYCLOCIM/Grid/M3d91x180x24.mat M3d grd;
  load /DFS-L/DATA/primeau/salali/CTL_OCIM_gas_exchange_OCT2023/DATA/M3d91x180x24.mat M3d grd;
  iocn = find(M3d==1); nocn = length(iocn);

% load the GLODAPv2 data
if fr == 1
  disp(sprintf('now loading %s from GLODAPv2.2022',vn)); 
  dpath = '/DFS-L/DATA/primeau/salali/CTL_OCIM_gas_exchange_OCT2023/DATA';
  %fpath = sprintf('%s/GLODAPv2/GLODAPv2.2022_Atlantic_Ocean.csv',dpath);
  fpath = sprintf('%s/GLODAPv2.2022_Atlantic_Ocean.csv',dpath);
  D{1} = importdata(fpath);
  %keyboard;
  fpath = sprintf('%s/GLODAPv2.2022_Pacific_Ocean.csv',dpath);
  D{2} = importdata(fpath);
  fpath = sprintf('%s/GLODAPv2.2022_Indian_Ocean.csv',dpath);
  D{3} = importdata(fpath);
  fpath = sprintf('%s/GLODAPv2.2022_Arctic_Ocean.csv',dpath);
  D{4} = importdata(fpath);
  Gdat = [D{1}.data;D{2}.data;D{3}.data;D{4}.data]; 
else
   disp(sprintf('now loading %s from GLODAPv2.2022',vn)); 
   dpath= '/DFS-L/DATA/primeau/salali/CTL_OCIM_gas_exchange_OCT2023/DATA';
   fpath=sprintf('%s/GLODAPv2.2022_Merged_Master_File.mat',dpath);
   load(fpath,'G2latitude','G2longitude','G2pressure','G2year','G2pressure','G2year','G2month','G2depth','G2c14','G2c14f');
end

  % geo variables
if fr==1
  ilat = find(strcmp(D{1}.textdata,'G2latitude'));
  ilon = find(strcmp(D{1}.textdata,'G2longitude'));
  ipres = find(strcmp(D{1}.textdata,'G2pressure'));
  iyear = find(strcmp(D{1}.textdata,'G2year'));
  imonth = find(strcmp(D{1}.textdata,'G2month'));
  idepth=find(strcmp(D{1}.textdata,'G2depth'));
  lat = Gdat(:,ilat);
  lon = Gdat(:,ilon);
  pres = Gdat(:,ipres);
  year = Gdat(:,iyear);
  month = Gdat(:,imonth);
  %depth = sw_dpth(pres,lat);
  depth=Gdat(:,idepth);   
  % fix negative longitudes
  ineg = find(lon<0);
  lon(ineg) = 360+lon(ineg);
else
  lat=G2latitude;
  lon=G2longitude;
  pressure=G2pressure;
  year=G2year;
  month=G2month;
  depth=G2depth;
  %fix negative longitudes
  ineg=find(lon<0);
  lon(ineg)=360+lon(ineg);
end  

  %% extract oceanic GLODAPv2 data (CFC11,CFC12,TCO2,...)
  if fr==1
  disp(sprintf('extracting %s from %i to %i',vn,parms.t0,parms.tf));
  icfc11 = find(strcmp(D{1}.textdata,vn));
  cfc11 = Gdat(:,icfc11);
  icfc11f = find(strcmp(D{1}.textdata,vnv));
  cfc11f = Gdat(:,icfc11f);
else 
  disp(sprintf('extracting %s from %i to %i',vn,parms.t0,parms.tf));
  cfc11=G2c14;
  cfc11f=G2c14f;
end 
  cfc11err = max(0.1*cfc11,0.05);
  idat = find(cfc11~=-9999 & depth>=0 & cfc11f==2 & year>=parms.t0 & year<=parms.tf);
  CFC11.cfc11 = cfc11(idat);
  CFC11.cfc11err = cfc11err(idat);
  % % remove outliers of GLODAPv2 data
  [~,TF] = rmoutliers(CFC11.cfc11,'percentiles',[0.1 99.9]);
  ibadV = find(TF==1);
  CFC11.cfc11(ibadV) = -9.999;
  % outlier removed
  CFC11.lon = lon(idat);
  CFC11.lat = lat(idat);
  CFC11.z = depth(idat);
  CFC11.yr = year(idat);
  CFC11.mo = month(idat);
  CFC11.id = CFC11.yr*12+CFC11.mo; 
  ID = unique(CFC11.id);  
  % bin to grd (year,month)
  for i = 1:length(ID)
      id = find(CFC11.id==ID(i));
      [mu,vard,n] = bin3d(CFC11.lon(id),CFC11.lat(id),CFC11.z(id),...
                          CFC11.cfc11(id),CFC11.cfc11err(id),...
                          grd.XT3d,grd.YT3d,grd.ZT3d);  
      % data.cfc11star(:,i) = mu(iocn);
      eval(sprintf('data.%sstar(:,i) = mu(iocn);',vnw));
      % data.varcfc11(:,i) = vard(iocn);
      eval(sprintf('data.var%s(:,i) = vard(iocn);',vnw));
      % data.ncfc11(:,i) = n(iocn);
      eval(sprintf('data.n%s(:,i) = n(iocn);',vnw));
  end
  eval(sprintf('data.%sid = ID;',vnw));
  % data.cfc11id = ID; % to retrieve the month and year, use
                           % mod(id,12)-> month, floor(id/12)
  eval(sprintf('data.%sz = CFC11.z;',vnw));
  % data.cfc11z = CFC11.z;

  %% create operator to compare the monthly GLODAPv2 data and model result
  % id1 = data.cfc11id; % id for the observed data
  id1 = eval(sprintf('data.%sid',vnw)); % id for the observed data
  II1 = ones(12,1); % months
  % II2 = ones(70,1); % years 
  II2 = ones(parms.tf-parms.t0+1,1); % years 
  % idd1 = kron(II1,[1940:2009]'*12) + kron([1:12]',II2);
  idd1 = kron(II1,[parms.t0:parms.tf]'*12) + kron([1:12]',II2);
  l1 = length(id1);
  l2 = length(idd1);
  H = sparse(l1,l2); 
  for i = 1:l1
      j = find(idd1 == id1(i));
      H(i,j) = 1;
  end  
  % tranfer CFC(n1,t1)(:) -> CFC(n1,t2)(:) 
  H1 = kron(H,speye(nocn)); 
  % step2: find grd
  cfc11star = eval(sprintf('data.%sstar',vnw));
  H2 = d0(cfc11star(:)~=-9.9990);
  ikeep = find(cfc11star(:)~=-9.9990);
  % H2 = d0(data.cfc11star(:)~=-9.9990);
  % ikeep = find(data.cfc11star(:)~=-9.9990);
  H2 = H2(:,ikeep)'; 
  % data.CFC12h1 = H1;
  % data.CFC12h2 = H2;
  eval(sprintf('data.%sh1 = H1;',upper(vnw)));
  eval(sprintf('data.%sh2 = H2;',upper(vnw)));

  % creat operator to compare annual GLODAPv2 data and model results
  % havg = @(v,dV) squeeze(nansum(v.*dV,[1,2]))./squeeze(nansum(v.*dV./v,[1,2]));

end
