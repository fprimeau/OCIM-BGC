function DOCclean = RemoveRef(par) ;

M3d  = par.M3d    ;
MSKS = par.MSKS   ;
DOC  = par.DOCobs ;

% Arctic Ocean
junk = M3d ;
junk(1:77,:,:) = 0 ;
iARC = find(junk(:)) ;
DOC(iARC) = DOC(iARC) - 47.1 ;

% Southern Ocean
junk = M3d ;
junk(21:end,:,:) = 0 ;
iSO = find(junk(:))  ;
DOC(iSO) = DOC(iSO) - 40.9 ;

% west N Atlantic
junk = par.MSKS.ATL ;
junk([1:45,78:end],:,:) = 0  ;
NATL  = junk ;
WNATL = NATL ;
WNATL(:,[1:20,wrapTo360(-40)/2:end],:) = 0 ;
ENATL = junk - WNATL ;

iWNAtl = find(WNATL(:)) ;
iENAtl = find(ENATL(:)) ;
DOC(iWNAtl) = DOC(iWNAtl) - 43.6 ;
DOC(iENAtl) = DOC(iENAtl) - 42.9 ;

% South Atlantic
junk = par.MSKS.ATL ;
junk([1:20,46:end],:,:) = 0 ;
iSAtl = find(junk(:)) ;
DOC(iSAtl) = DOC(iSAtl) - 40.5 ;

% North Pacific
junk = par.MSKS.PAC   ;
junk(1:45,:,:) = 0   ;
iNPac = find(junk(:)) ;
DOC(iNPac) = DOC(iNPac) - 38.5 ;

% South Pacific
junk = par.MSKS.PAC   ;
junk([1:20,46:end],:,:) = 0 ;
iSPac = find(junk(:)) ;
DOC(iSPac) = DOC(iSPac) - 38.7 ;

% North Indian
junk = par.MSKS.IND   ;
junk(1:45,:,:) = 0   ;
iNInd = find(junk(:)) ;
DOC(iNInd) = DOC(iNInd) -41.7 ;

% South Indian
junk = par.MSKS.IND   ;
junk(46:end,:,:) = 0  ;
iSInd = find(junk(:)) ;
DOC(iSInd) = DOC(iSInd) - 41.2 ;

% Med.
Med = par.MSKS.MED ;
iMed = find(Med(:)) ;
DOC(iMed) = nan ;
DOCclean  = DOC ;
