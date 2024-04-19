%
% Observed DIC13 and DIC14 from GLODAPv2.2023 to OCIM2 24layer grid
%
addpath('../src')
addpath('../../Observation_rawdata/GLODAP')
addpath('../../DATA/BGC_2023Nature/')
load GLODAPv2_2023_Merged_Master_File.mat 
%https://glodap.info/index.php/merged-and-adjusted-data-product-v2-2023/
%Download GLODAPv2.2023 data from above link

%
load OCIM2_CTL_He.mat
M3d = output.M3d  ;
grd = output.grid ; 
x = wrapTo360(G2longitude) ;
y = G2latitude ;
z = G2depth    ;

%excluding NaN from GLODAP data and use bin3d to bin to the OCIM grid.

% c13f: WOCE flag name
% c13qc: Secondary flag name
% c14err: Counting error
% delC13:   [permil]
% deltaC14: [permil]

% See following Reference--> The annual update GLODAPv2.2023: the global interior ocean biogeochemical data product, 2024, Earth System Science Data.

ic13data = find(~isnan(G2c13));
c13lon   = x(ic13data);
c13lat   = y(ic13data);
c13dep   = z(ic13data);
c13raw   = G2c13(ic13data);
c13err   = c13raw*0.01;

ic14data = find(~isnan(G2c14));
c14lon   = x(ic14data);
c14lat   = y(ic14data);
c14dep   = z(ic14data);
c14raw   = G2c14(ic14data);
c14err   = c14raw*0.01;

%
[c13raw,c13err]     = bin3d(c13lon,c13lat,c13dep,c13raw,c13err,grd.XT3d,grd.YT3d,grd.ZT3d);
[c14raw,c14err]     = bin3d(c14lon,c14lat,c14dep,c14raw,c14err,grd.XT3d,grd.YT3d,grd.ZT3d);

%
c13raw(c13raw(:) ==-9.999)=NaN;
c14raw(c14raw(:) ==-9.999)=NaN;
%
c13raw(M3d(:)==0)=NaN;
c14raw(M3d(:)==0)=NaN;


%DATA file deirectory에 저장하기
fileName = 'GLODAPv2_DICisotopes_91x180x24.mat'
directory = '../../DATA/BGC_2023Nature'
filePath = fullfile(directory, fileName);
save(filePath, 'c13raw', 'c14raw');
