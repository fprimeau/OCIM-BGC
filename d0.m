function X = d0(v) 
X = spdiags(v(:),0,length(v(:)),length(v(:)));
