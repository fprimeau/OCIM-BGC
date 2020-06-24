using Plots, MAT, MATLAB
matfile1 = matopen("/DFS-L/DATA/primeau/weilewang/DATA/transport_v4.mat");
matfile2 = matopen("/DFS-L/DATA/primeau/weilewang/DATA/tempobs_90x180x24.mat");

Temp = read(matfile2,"tempobs");
@mput Temp
eval_string("aveT = nanmean(Temp(:,:,1:8),3)");
@mget aveT;

grd = read(matfile1,"grid");
x = grd["xt"][:];
y = grd["yt"][:];

αₚ = -9.99e-03;
βₚ = 1.09;

αₒ = 7.18e-03;
βₒ = 1.02;

bₚ = αₚ*aveT .+ βₚ;
bₒ = αₒ*aveT .+ βₒ;

title1 = "b distribution for P";
title2 = "b distribution for C";
plt1 = contourf(x,y,bₚ,title = title1);
plt2 = contourf(x,y,bₒ, title = title2);
plot(plt1,plt2,layout = (2,1));
png("b_distribution.pdf")
