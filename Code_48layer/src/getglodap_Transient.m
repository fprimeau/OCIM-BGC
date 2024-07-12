function data = getglodap_Transient(t0,tf)
    % function data = getglodap_c14(t0,tf)
    %     data: required, struct  
    %       vn: required, variable name for extract
    %       t0: optional, starting year for data extraction
    %       tf: optional, ending year for data extraction
    % extract variable from GLODAPv2.2021 data
    % interpolate data to OCIM grid
    % build an operator to map OCIM data to where
    % GLODAP data available
    
    load /DFS-L/DATA/primeau/hojons1/OCIM-BGC/DATA/BGC_48layer/OCIM2_CTL_He_48layer.mat output % call the M3d and output
    M3d = output.M3d  ;
    grd = output.grid ;
    iwet = find(M3d==1); 
    nwet = length(iwet);
    
    % load the GLODAPv2 data
    disp('now loading C13 and C14 data from GLODAPv2.2023'); 
    
    dpath= '/DFS-L/DATA/primeau/hojons1/OCIM-BGC/Observation_rawdata/GLODAP/';
    fpath=sprintf('%s/GLODAPv2_2023_Merged_Master_File.mat',dpath);
    load(fpath,'G2latitude','G2longitude','G2pressure','G2year','G2month','G2depth','G2tco2','G2tco2f','G2c13','G2c13f','G2c14','G2c14err','G2c14f');
    
    lat=G2latitude;
    lon=G2longitude;
    pressure=G2pressure;
    year=G2year;
    month=G2month;
    depth=G2depth;
    %fix negative longitudes
    ineg=find(lon<0);
    lon(ineg)=360+lon(ineg);
    
    
    disp(sprintf('extracting C13(permil) and C14(permil) from %i to %i',t0,tf));
    dic  = G2tco2   ;
    dicf = G2tco2f  ;
    c13  = G2c13     ;
    c13f = G2c13f    ;
    c14  = G2c14     ;
    c14f = G2c14f    ;
    %o2   = G2oxygen  ;
    %o2f  = G2oxygenf ;
    
    dicerr = max(0.005*dic, 10);
    c13err = max(0.1*c13,0.05); % NOT CORRECT
    c14err = G2c14err ;
    %o2err  = max(0.005*c13,1); % NOT CORRECT
    %
    idatdic = find( ( dic ~= -9999 ) & ( depth >= 0 ) & ( dicf==2 ) & ( year >= t0 ) & ( year <= tf) );
    dic = dic(idatdic);
    dicerr = dicerr(idatdic);
    %
    idatc13 = find( ( c13 ~= -9999 ) & ( depth >= 0 ) & ( c13f==2 ) & ( year >= t0 ) & ( year <= tf) );
    c13 = c13(idatc13);
    c13err = c13err(idatc13);
    %    
    idatc14 = find( ( c14 ~= -9999 ) & ( depth >= 0 ) & ( c14f==2 ) & ( year >= t0 ) & ( year <= tf) );
    c14 = c14(idatc14);
    c14err = c14err(idatc14);
    %    
    %idato2 = find( ( o2 ~= -9999 ) & ( depth >= 0 ) & ( o2f==2 ) & ( year >= t0 ) & ( year <= tf) );
    %o2 = o2(idato2);
    %o2err = o2err(idato2);
    
    %
    londic = lon(idatdic);
    latdic = lat(idatdic);
    zdic = depth(idatdic);
    yrdic = year(idatdic);
    modic = month(idatdic);
    modic = floor((modic-1)/2) + 1 ; % for the 2-month time-step.
    idxdic = yrdic * 6 + modic;  % you can adjust depending on the time interval
    IDdic = unique(idxdic); 
    %
    lonc13 = lon(idatc13);
    latc13 = lat(idatc13);
    zc13 = depth(idatc13);
    yrc13 = year(idatc13);
    moc13 = month(idatc13);
    moc13 = floor((moc13-1)/2) + 1 ;
    idxc13 = yrc13 * 6 + moc13;  % you can adjust depending on the time interval
    IDc13 = unique(idxc13); 
    %
    lonc14 = lon(idatc14);
    latc14 = lat(idatc14);
    zc14 = depth(idatc14);
    yrc14 = year(idatc14);
    moc14 = month(idatc14);
    moc14 = floor((moc14-1)/2) + 1 ;
    idxc14 = yrc14 * 6 + moc14;  % you can adjust depending on the time interval
    IDc14 = unique(idxc14); 
    %
    %lono2 = lon(idato2);
    %lato2 = lat(idato2);
    %zo2 = depth(idato2);
    %yro2 = year(idato2);
    %moo2 = month(idato2);
    %idxo2 = yro2 * 12 + moo2; 
    %IDo2 = unique(idxo2); 

    
    % bin to grd (year,month)
    %----------for dic--------------------------
    for i = 1:length(IDdic)
        iddic = find( idxdic == IDdic(i) );
        
        [mu,vard,n] = bin3d(londic(iddic), latdic(iddic), zdic(iddic), dic(iddic), dicerr(iddic), grd.XT3d, grd.YT3d, grd.ZT3d);  
        
        data.dic(:,i)    = mu(iwet);
        data.vardic(:,i) = vard(iwet);
        data.ndic(:,i)   = n(iwet);
    end
    data.dicid = IDdic;
    data.dicz  = zdic;
    %----------for C13--------------------------
    for i = 1:length(IDc13)
        idc13 = find( idxc13 == IDc13(i) );
        
        [mu,vard,n] = bin3d(lonc13(idc13), latc13(idc13), zc13(idc13), c13(idc13), c13err(idc13), grd.XT3d, grd.YT3d, grd.ZT3d);  
        
        data.c13(:,i)    = mu(iwet);
        data.varc13(:,i) = vard(iwet);
        data.nc13(:,i)   = n(iwet);
    end
    data.c13id = IDc13;
    data.c13z  = zc13;
    %------------for C14----------------------
    for i = 1:length(IDc14)
        idc14 = find( idxc14 == IDc14(i) );
        
        [mu,vard,n] = bin3d(lonc14(idc14), latc14(idc14), zc14(idc14), c14(idc14), c14err(idc14), grd.XT3d, grd.YT3d, grd.ZT3d);  
        
        data.c14(:,i)    = mu(iwet);
        data.varc14(:,i) = vard(iwet);
        data.nc14(:,i)   = n(iwet);
    end
    data.c14id = IDc14;
    data.c14z  = zc14;
    %------------O2----------------------
    %for i = 1:length(IDo2)
    %    ido2 = find( idxo2 == IDo2(i) );
         
    %    [mu,vard,n] = bin3d(lono2(ido2), lato2(ido2), zo2(ido2), o2(ido2), o2err(ido2), grd.XT3d, grd.YT3d, grd.ZT3d);  
            
    %    data.o2(:,i)    = mu(iwet);
    %    data.varo2(:,i) = vard(iwet);
    %    data.no2(:,i)   = n(iwet);
    %end
    %data.o2id = IDo2;
    %data.o2z  = zo2;
 
    %% create operator to compare the monthly GLODAPv2 data and model result
    %----------for DIC-------------------
    id1_dic = data.dicid;
    idd1 = [];
    for i = t0:tf
        for j = 1:6
            idd1 = [idd1; i*6 + j];
        end
    end
    l1_dic = length(id1_dic);
    l2 = length(idd1);
    H_dic = sparse(l1_dic,l2); 
    for i = 1:l1_dic
        j = find(idd1 == id1_dic(i));
        H_dic(i,j) = 1;
    end  
    H1_dic = kron(H_dic,speye(nwet));
    %----------for C13-------------------
    id1_c13 = data.c13id;
    idd1 = [];
    for i = t0:tf
        for j = 1:6
            idd1 = [idd1; i*6 + j];
        end
    end
    l1_c13 = length(id1_c13);
    l2 = length(idd1);
    H_c13 = sparse(l1_c13,l2); 
    for i = 1:l1_c13
        j = find(idd1 == id1_c13(i));
        H_c13(i,j) = 1;
    end  
    H1_c13 = kron(H_c13,speye(nwet)); 
    %----------for C14-------------------
    id1_c14 = data.c14id;
    idd1 = [];
    for i = t0:tf
        for j = 1:6
            idd1 = [idd1; i*6 + j];
        end
    end
    l1_c14 = length(id1_c14);
    l2 = length(idd1);
    H_c14 = sparse(l1_c14,l2); 
    for i = 1:l1_c14
        j = find(idd1 == id1_c14(i));
        H_c14(i,j) = 1;
    end  
    H1_c14 = kron(H_c14,speye(nwet));
    %----------for O2-------------------
    %id1_o2 = data.o2id;
    %II1 = ones(12,1); % months
    %II2 = ones(tf-t0+1,1); % years 
    %idd1 = kron(II1,[t0:tf]'*12) + kron([1:12]',II2);
    %l1_o2 = length(id1_o2);
    %l2 = length(idd1);
    %H_o2 = sparse(l1_o2,l2); 
    %for i = 1:l1_o2
    %    j = find(idd1 == id1_o2(i));
    %    H_o2(i,j) = 1;
    %end  
    %H1_o2 = kron(H_o2,speye(nwet));  


    % step2: find grd
    dic = data.dic;
    d0 = @(x) spdiags(x(:),[0],length(x(:)),length(x(:)));
    H2_dic = d0(dic(:)~=-9.999);
    ikeep_dic = find(dic(:)~=-9.999);
    H2_dic = H2_dic(:,ikeep_dic)'; 
    %----------for C13--------------
    c13 = data.c13;
    d0 = @(x) spdiags(x(:),[0],length(x(:)),length(x(:)));
    H2_c13 = d0(c13(:)~=-9.999);
    ikeep_c13 = find(c13(:)~=-9.999);
    H2_c13 = H2_c13(:,ikeep_c13)'; 
    %----------for C14--------------
    c14 = data.c14;
    d0 = @(x) spdiags(x(:),[0],length(x(:)),length(x(:)));
    H2_c14 = d0(c14(:)~=-9.999);
    ikeep_c14 = find(c14(:)~=-9.999);
    H2_c14 = H2_c14(:,ikeep_c14)'; 
    %----------for O2--------------
    %o2 = data.o2;
    %d0 = @(x) spdiags(x(:),[0],length(x(:)),length(x(:)));
    %H2_o2 = d0(o2(:)~=-9999);
    %ikeep_o2 = find(o2(:)~=-9999);
    %H2_o2 = H2_o2(:,ikeep_o2)'; 
   
    %
    data.dich1   = H1_dic;
    data.dich2   = H2_dic;
    data.dic13h1 = H1_c13;
    data.dic13h2 = H2_c13;
    data.dic14h1 = H1_c14;
    data.dic14h2 = H2_c14;
    %data.o2h1    = H1_o2 ;
    %data.o2h2    = H2_o2 ;
    
end
