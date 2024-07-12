function cleanD = rmOutliers(obs, M3d) ;
    iwet = find(M3d(:)) ;
    vD = obs(iwet) ;
    ineg = find(vD<0) ;
    vD(ineg) = nan    ;
    [~,TF] = rmoutliers(vD,'percentiles',[0.1 99.9]);
    ibadD = find(TF==1) ;
    vD(ibadD) = nan ;
    cleanD = M3d+nan;
    cleanD(iwet) = vD ;
end

