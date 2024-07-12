function [S, I, rl, rs, rm] = regression(data1, data2)
%
% calculate the slope(S), intercept(I), regression line(rl) r-square(rs), rmse between data1 (X-axis) and data2 (Y-axis) (linear regression)
% if data has an nan values, nana values should be removed first
% inan = isnan(dic13raw) | isnan(delDIC13) ;
% e.g., data1: Observational Results; data2: Model Results

fitResult = polyfit(data1, data2, 1) ;
S = fitResult(1) ;
I = fitResult(2) ;
%
rl   = polyval(fitResult, data1) ;
yMean = mean(data2) ;
ssTot = sum((data2 - yMean).^2) ;
ssRes = sum((data2 - rl).^2) ;
%
rs = 1-(ssRes/ssTot) ;
rm = rmse(data2, data1) ;

end
