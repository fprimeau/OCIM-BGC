using Plots, MAT, MATLAB
mf1 = matopen("/DFS-L/DATA/primeau/weilewang/DATA/po4obs_90x180x24.mat");
mf2 = matopen("/DFS-L/DATA/primeau/weilewang/DATA/transport_v4.mat");

po4obs = read(mf1, "po4obs");
grd = read(mf2, "grid");
x = grd["xt"][:];
y = grd["yt"][:];

DIP = po4obs[:,:,1];
cc = 1.36e-03 ;
dd = 6.38e-03 ;

P2C = cc*DIP .+ dd;
C2P = 1 ./ P2C;

title = "C:P ratio";
# pyplot()
plt = contourf(x, y, C2P,title = title);
png("C2P.png")
