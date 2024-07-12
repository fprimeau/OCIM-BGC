function data = getOC_Transient(t0,tf)
    % function data = getglodap_c14(t0,tf)
    %     data: required, struct  
    %       vn: required, variable name for extract
    %       t0: optional, starting year for data extraction
    %       tf: optional, ending year for data extraction
    % extract variable from GLODAPv2.2021 data
    % interpolate data to OCIM grid
    % build an operator to map OCIM data to where
    % GLODAP data available
    
    addpath('../../Observation_rawdata/')
    load /DFS-L/DATA/primeau/hojons1/OCIM-BGC/DATA/BGC_48layer/OCIM2_CTL_He_48layer.mat output % call the M3d and output
    M3d = output.M3d  ;
    grd = output.grid ;
    

    iwet = find(M3d==1); nwet = length(iwet);
    % load the DO13C and DO14C data
    disp('now loading DOC13 and DOC14 data from Hansell et al'); 
    dpath= '/DFS-L/DATA/primeau/hojons1/OCIM-BGC/Observation_rawdata/';
    fpath=sprintf('%s/DOM_Merged_Hansell_2022_Cisotopes.xlsx',dpath);
    dataTable = readtable(fpath, 'Sheet', 'MergedBasins_V5', 'VariableNamingRule', 'preserve') ;
    % remove wrong data format (I also check there is no C13 and C14 data in this weird data format)
    dataTable(1,:) = [];
    date  = table2array(dataTable(:,1)) ; %YYYYMMDD
    datestr = string(date) ;
    nvdata = strlength(datestr) ~= 8 ;
    invdata = find(nvdata == 1);
    dataTable(invdata, :) = [];

    %
    londoc   = table2array(dataTable(:,3)); %o-180ow to 180oE. Thus, we need wrapTo360
    latdoc   = table2array(dataTable(:,2));
    depthdoc = table2array(dataTable(:,5)); %In rawdata: CTD pressure (dbar)
    doc13 = table2array(dataTable(:,25)); %permil
    doc14 = table2array(dataTable(:,27)); %permil
    date  = table2array(dataTable(:,1)) ; %YYYYMMDD
    datestr = num2str(date, '%08d') ;
    yeardoc    = str2num(datestr(:,1:4)) ;
    monthdoc   = str2num(datestr(:,5:6)) ;
    %fix negative longitudes
    ineg=find(londoc<0);
    londoc(ineg)=360+londoc(ineg);

    % load additional isotope data from Hojong's compiliation-----
    % References are in the excel data file.
    disp('now loading DOC13 and DOC14 data from mannually compiled excel files');
    dataTable_compil = readtable('DOCisotope_compilation_Hojong.xlsx', 'Sheet', 'Sheet1', 'VariableNamingRule', 'preserve');
    dataTable_compil(1,:) = [];
    lon_compil      = table2array(dataTable_compil(:,5)); %0oE to 360oE. 
    lat_compil      = table2array(dataTable_compil(:,4));
    dep_compil      = table2array(dataTable_compil(:,6)); %
    doc13_compil    = table2array(dataTable_compil(:,7)); %permil
    doc14_compil    = table2array(dataTable_compil(:,8)); %permil
    yeardoc_compil  = table2array(dataTable_compil(:,2));
    monthdoc_compil = table2array(dataTable_compil(:,3));

    %-------------combined DOC13 and DOC14 data from two diffent references-----------------
    latdoc = [latdoc; lat_compil] ;
    londoc = [londoc; lon_compil] ;
    depthdoc = [depthdoc; dep_compil] ;
    yeardoc = [yeardoc; yeardoc_compil] ;
    monthdoc = [monthdoc; monthdoc_compil] ;
    doc13 = [doc13; doc13_compil] ;
    doc14 = [doc14; doc14_compil] ;

    % load the po13C data
    disp('now loading poc13 data from Verwega et al. (2021)');
    % data from https://doi.pangaea.de/10.1594/PANGAEA.946915?format=html#download
    % Detail information can be referred in following paper
    % Description of a global marine particulate organic carbon-13 isotope data set, Verwega et al. (2021), Earth System Science Data
    poc13    = ncread('poc13_all_year1966-2019_woa09grid.nc','POC13') ; %del13C-POC permil, missing value:-9999
    timepoc  = ncread('poc13_all_year1966-2019_woa09grid.nc','TIME')  ;
    depthpoc = ncread('poc13_all_year1966-2019_woa09grid.nc','DEPTH') ;
    latpoc   = ncread('poc13_all_year1966-2019_woa09grid.nc','LAT')   ;
    lonpoc   = ncread('poc13_all_year1966-2019_woa09grid.nc','LON')   ;
    poc13  =  permute(poc13,   [2, 1, 3, 4 ]); 
    [latpoc, lonpoc, depthpoc, timepoc]    = ndgrid(latpoc, lonpoc, depthpoc, timepoc) ;
    %
    base_date = datetime(1966, 6, 16, 12, 0, 0);
    date = base_date + days(timepoc- 1);
    yearpoc =  year(date);
    monthpoc = month(date);

   
    %-------------find data point-------------------
    doc13err  = doc13*0.02; % NOT CORRECT
    doc14err =  doc14*0.02; % NOT CORRECT
    poc13err =  poc13*0.02; % NOT CORRECT
    %
    idatdoc13 = find( ~isnan(doc13) & ( doc13 ~= -999 ) & ( depthdoc >= 0 ) & ( yeardoc >= t0 ) & ( yeardoc <= tf) );
    doc13    = doc13(idatdoc13);
    doc13err = doc13err(idatdoc13);
    %    
    idatdoc14 = find( ~isnan(doc14) & ( doc14 ~= -999 ) & ( depthdoc >= 0 ) & ( yeardoc >= t0 ) & ( yeardoc <= tf) );
    doc14     = doc14(idatdoc14);
    doc14err  = doc14err(idatdoc14);
    %    
    idatpoc13 = find( ~isnan(poc13) & ( depthpoc >= 0 ) & ( yearpoc >= t0 ) & ( yearpoc <= tf) );
    poc13     = poc13(idatpoc13);
    poc13err  = poc13err(idatpoc13);
    
    %
    londoc13 = londoc(idatdoc13);
    latdoc13 = latdoc(idatdoc13);
    zdoc13   = depthdoc(idatdoc13);
    yrdoc13  = yeardoc(idatdoc13);
    modoc13  = monthdoc(idatdoc13);
    modoc13 = floor((modoc13-1)/2) + 1 ; % for the 2-month time-step. It makes Jan-Feb to 1, Mar-Apr to 2, ... Nov-Dec to 6.
    idxdoc13 = yrdoc13 * 6 + modoc13;    % you can adjust depending on the time interval
    IDdoc13  = unique(idxdoc13); 
    %
    londoc14 = londoc(idatdoc14);
    latdoc14 = latdoc(idatdoc14);
    zdoc14   = depthdoc(idatdoc14);
    yrdoc14  = yeardoc(idatdoc14);
    modoc14  = monthdoc(idatdoc14);
    modoc14 = floor((modoc14-1)/2) + 1 ;
    idxdoc14 = yrdoc14 * 6 + modoc14;  % you can adjust depending on the time interval
    IDdoc14  = unique(idxdoc14); 
    %
    lonpoc13 = lonpoc(idatpoc13);
    latpoc13 = latpoc(idatpoc13);
    zpoc13   = depthpoc(idatpoc13);
    yrpoc13  = yearpoc(idatpoc13);
    mopoc13  = monthpoc(idatpoc13);
    mopoc13 = floor((mopoc13-1)/2) + 1 ;
    idxpoc13 = yrpoc13 * 6 + mopoc13;  % you can adjust depending on the time interval
    IDpoc13  = unique(idxpoc13); 

    
    % bin to grd (year,month)
    %----------for DOC13--------------------------
    for i = 1:length(IDdoc13)
        iddoc13 = find( idxdoc13 == IDdoc13(i) );
        
        [mu,vard,n] = bin3d(londoc13(iddoc13), latdoc13(iddoc13), zdoc13(iddoc13), doc13(iddoc13), doc13err(iddoc13), grd.XT3d, grd.YT3d, grd.ZT3d);  
        
        data.doc13(:,i)    = mu(iwet);
        data.vardoc13(:,i) = vard(iwet);
        data.ndoc13(:,i)   = n(iwet);
    end
    data.doc13id = IDdoc13;
    data.doc13z  = zdoc13;
    %------------for DOC14----------------------
    for i = 1:length(IDdoc14)
        iddoc14 = find( idxdoc14 == IDdoc14(i) );
        
        [mu,vard,n] = bin3d(londoc14(iddoc14), latdoc14(iddoc14), zdoc14(iddoc14), doc14(iddoc14), doc14err(iddoc14), grd.XT3d, grd.YT3d, grd.ZT3d);  
        
        data.doc14(:,i)    = mu(iwet);
        data.vardoc14(:,i) = vard(iwet);
        data.ndoc14(:,i)   = n(iwet);
    end
    data.doc14id = IDdoc14;
    data.doc14z  = zdoc14;
   %----------for POC13--------------------------
   for i = 1:length(IDpoc13)
    idpoc13 = find( idxpoc13 == IDpoc13(i) );
    
    [mu,vard,n] = bin3d(lonpoc13(idpoc13), latpoc13(idpoc13), zpoc13(idpoc13), poc13(idpoc13), poc13err(idpoc13), grd.XT3d, grd.YT3d, grd.ZT3d);  
    
    data.poc13(:,i)    = mu(iwet);
    data.varpoc13(:,i) = vard(iwet);
    data.npoc13(:,i)   = n(iwet);
    end
    data.poc13id = IDpoc13;
    data.poc13z  = zpoc13;  
 
    %% create operator to compare the monthly GLODAPv2 data and model result
    %----------for DOC13-------------------
    id1_doc13 = data.doc13id;
    idd1 = [];
    for i = t0:tf
        for j = 1:6
            idd1 = [idd1; i*6 + j];
        end
    end
    l1_doc13 = length(id1_doc13);
    l2 = length(idd1);
    H_doc13 = sparse(l1_doc13,l2); 
    for i = 1:l1_doc13
        j = find(idd1 == id1_doc13(i));
        H_doc13(i,j) = 1;
    end  
    H1_doc13 = kron(H_doc13,speye(nwet)); 
    %----------for DOC14-------------------
    id1_doc14 = data.doc14id;
    idd1 = [];
    for i = t0:tf
        for j = 1:6
            idd1 = [idd1; i*6 + j];
        end
    end
    l1_doc14 = length(id1_doc14);
    l2 = length(idd1);
    H_doc14 = sparse(l1_doc14,l2); 
    for i = 1:l1_doc14
        j = find(idd1 == id1_doc14(i));
        H_doc14(i,j) = 1;
    end  
    H1_doc14 = kron(H_doc14,speye(nwet));
    %----------for POC13-------------------
    id1_poc13 = data.poc13id;
    idd1 = [];
    for i = t0:tf
        for j = 1:6
            idd1 = [idd1; i*6 + j];
        end
    end
    l1_poc13 = length(id1_poc13);
    l2 = length(idd1);
    H_poc13 = sparse(l1_poc13,l2); 
    for i = 1:l1_poc13
        j = find(idd1 == id1_poc13(i));
        H_poc13(i,j) = 1;
    end  
    H1_poc13 = kron(H_poc13,speye(nwet));  


    % step2: find grd
    %----------for DOC13--------------
    doc13 = data.doc13;
    d0 = @(x) spdiags(x(:),[0],length(x(:)),length(x(:)));
    H2_doc13 = d0(doc13(:)~=-9.999);
    ikeep_doc13 = find(doc13(:)~=-9.999);
    H2_doc13 = H2_doc13(:,ikeep_doc13)'; 
    %----------for DOC14--------------
    doc14 = data.doc14;
    d0 = @(x) spdiags(x(:),[0],length(x(:)),length(x(:)));
    H2_doc14 = d0(doc14(:)~=-9.999);
    ikeep_doc14 = find(doc14(:)~=-9.999);
    H2_doc14 = H2_doc14(:,ikeep_doc14)'; 
    %----------for O2--------------
    poc13 = data.poc13;
    d0 = @(x) spdiags(x(:),[0],length(x(:)),length(x(:)));
    H2_poc13 = d0(poc13(:)~=-9.999);
    ikeep_poc13 = find(poc13(:)~=-9.999);
    H2_poc13 = H2_poc13(:,ikeep_poc13)'; 
   
    %
    data.doc13h1 = H1_doc13;
    data.doc13h2 = H2_doc13;
    data.doc14h1 = H1_doc14;
    data.doc14h2 = H2_doc14;
    data.poc13h1  = H1_poc13 ;
    data.poc13h2  = H2_poc13 ;
    
end
