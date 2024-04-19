%
%
addpath('../src')
addpath('../../Observation_rawdata/')
addpath('../../DATA/BGC_2023Nature/')

%
dataTable = readtable('DOM_Merged_Hansell_2022_Cisotopes.xlsx', 'Sheet', 'MergedBasins_V5', 'VariableNamingRule', 'preserve');
%https://www.ncei.noaa.gov/access/ocean-carbon-acidification-data-system/oceans/ndp_109/ndp109.html
%Download from above link and edit for DO13C and DO14C (raw file is also downloaded)

dataTable(1,:) = [];
longtitude   = table2array(dataTable(:,3)); %o-180ow to 180oE. Thus, we need wrapTo360
latitude     = table2array(dataTable(:,2));
depth        = table2array(dataTable(:,5)); %In rawdata: CTD pressure (dbar) ----> depth로 바꾸긴 해야함.
DOC13        = table2array(dataTable(:,25)); %permil
DOC14        = table2array(dataTable(:,27)); %permil

load OCIM2_CTL_He.mat
M3d = output.M3d  ;
grd = output.grid ; % grid >> XT3d, YT3d, ZT3d. 3d array (91 by 180 by 24)

x = wrapTo360(longtitude) ;
y = latitude ;
z = depth    ;

doc(doc <= 0) = NaN;
idocdata = find(~isnan(doc));
doclon   = x(idocdata);
doclat   = y(idocdata);
docdep   = z(idocdata);
docraw   = doc(idocdata);
docerr   = docraw*0.05;      % No reports on analytical uncertainty in REF. Shoule be checked.


[docraw,docerr]     = bin3d(doclon,doclat,docdep,docraw,docerr,grd.XT3d,grd.YT3d,grd.ZT3d);

docraw(docraw(:)==-9.999)=NaN;
docraw(M3d(:)==0)=NaN;

%DATA file deirectory에 저장하기
fileName = 'docraw_91x180x48.mat'
directory = '../../DATA/BGC_48layer'
filePath = fullfile(directory, fileName);
save(filePath, 'docraw');