function [D14c,d13c,DI12C,DI13C,DI14C,ddo13c,dpo13c,Ddo14c,Dpo14c]=load_results(path,fxkw)
    nwet   = 200160; % number of wet points in ocean
    nfile = 172; % number of files in the model output directory
    nstep = 6; % number of time steps per file
    tstep = nfile * nstep ;
    
    D14c = zeros(nwet, tstep); 
    DI14C= zeros(nwet, tstep); 
    DIC  = zeros(nwet, tstep); 
    d13c = zeros(nwet, tstep); %200160*1032
    ddo13c = zeros(nwet, tstep); %200160*1032
    dpo13c = zeros(nwet, tstep); %200160*1032
    Ddo14c = zeros(nwet, tstep); %200160*1032
    Dpo14c = zeros(nwet, tstep); %200160*1032
    
    
    for i=1:172 %change length to number of files in the model output directory
            fname=sprintf('Update2_Transient_%.4fxkw_fras=1.00_frpho=0.50_fc14=2.00_%i.mat',fxkw,i);
            s=[path '/' fname];
            
            load(s)
    
            % solve for Delta C14
            R14oxa=1.176*1e-12;
            R13pdb=1.12372*1e-2;
            R13 = @(c13,c12) c13./c12;
            R14 = @(c14,c12) c14./c12;
            delta13c = @(c13,c12) (( R13(c13,c12) ./  R13pdb ) - 1)  * 1e3;
            delta14c = @(c14,c12) (( R14(c14,c12) ./  R14oxa ) - 1)  * 1e3;
            Delta14c = @(c14,c13,c12) (( R14(c14,c12) ./ R14oxa ) .* ( 0.975 ./ ( 1 + delta13c(c13,c12) / 1e3 ) ).^2 -1 ) * 1e3;
            %          
            DIC12 =Xout.C12(0*nwet+1:1*nwet,:) - Xout.C13(0*nwet+1:1*nwet,:) ;
            POC12 =Xout.C12(1*nwet+1:2*nwet,:) - Xout.C13(1*nwet+1:2*nwet,:) ; 
            DOC12 =Xout.C12(2*nwet+1:3*nwet,:) - Xout.C13(2*nwet+1:3*nwet,:) ; 
            DOC12l=Xout.C12(5*nwet+1:6*nwet,:) - Xout.C13(4*nwet+1:5*nwet,:) ; 
            DOC12r=Xout.C12(6*nwet+1:7*nwet,:) - Xout.C13(5*nwet+1:6*nwet,:) ; 
            DOC12t = DOC12 + DOC12l + DOC12r ;       
            %          
            DIC13 =Xout.C13(1:nwet,:);
            POC13 =Xout.C13(1*nwet+1:2*nwet,:);
            DOC13 =Xout.C13(2*nwet+1:3*nwet,:);
            DOC13l=Xout.C13(4*nwet+1:5*nwet,:);
            DOC13r=Xout.C13(5*nwet+1:6*nwet,:);
            DOC13t = DOC13 + DOC13l + DOC13r ;
            %          
            DIC14 =Xout.C14(1:nwet,:);
            POC14 =Xout.C14(1*nwet+1:2*nwet,:);
            DOC14 =Xout.C14(2*nwet+1:3*nwet,:);
            DOC14l=Xout.C14(4*nwet+1:5*nwet,:);
            DOC14r=Xout.C14(5*nwet+1:6*nwet,:);
            DOC14t = DOC14 + DOC14l + DOC14r ;
            %
            tmp1= Delta14c(DIC14,DIC13,DIC12); %200160 x 6
            tmp2= delta13c(DIC13,DIC12);
            tmp3= DIC12 ;
            tmp4= DIC13 ;
            tmp5= DIC14 ;
            tmp6= delta13c(DOC13t,DOC12t) ; %delta DOC13
            tmp7= delta13c(POC13,POC12) ; %delta POC13
            tmp8= Delta14c(DOC14t,DOC13t,DOC12t) ; %Delta DOC14
            tmp9= Delta14c(POC14,POC13,POC12) ; %Delta POC14
           
            %
            istart =  (i-1) * nstep + 1 ;
            iend   = istart + nstep - 1 ;
            D14c (:, istart:iend)  = tmp1;
            d13c (:, istart:iend)  = tmp2;
            DI12C(:, istart:iend)  = tmp3;
            DI13C(:, istart:iend)  = tmp4;
            DI14C(:, istart:iend)  = tmp5;
            ddo13c(:, istart:iend) = tmp6;
            dpo13c(:, istart:iend) = tmp7;
            Ddo14c(:, istart:iend) = tmp8;
            Dpo14c(:, istart:iend) = tmp9;
    
    end
    end
    
    
    
    
    
    