%% Voronoi disk
addpath('../Geometry');
cx = 0;
cy = 0;
r = 1;

N = 2500;
hs = haltonset(2);
x = (cx-r) + 2*r*hs(1:N,1);
y = (cy-r) + 2*r*hs(1:N,2);
clipping_region = nsidedpoly(round(sqrt(N)*pi),'Center',[cx,cy],'Radius',r);
p = polysoup_from_voronoi_uniform(x,y,clipping_region,2);
[v,e,c] = polymesh_from_polysoup(p);
polymesh_save('voronoi_disk_2500',v,e,c);

N = 10000;
hs = haltonset(2);
x = (cx-r) + 2*r*hs(1:N,1);
y = (cy-r) + 2*r*hs(1:N,2);
clipping_region = nsidedpoly(round(sqrt(N)*pi),'Center',[cx,cy],'Radius',r);
p = polysoup_from_voronoi_uniform(x,y,clipping_region,2);
[v,e,c] = polymesh_from_polysoup(p);
polymesh_save('voronoi_disk_10000',v,e,c);

N = 40000;
hs = haltonset(2);
x = (cx-r) + 2*r*hs(1:N,1);
y = (cy-r) + 2*r*hs(1:N,2);
clipping_region = nsidedpoly(round(sqrt(N)*pi),'Center',[cx,cy],'Radius',r);
p = polysoup_from_voronoi_uniform(x,y,clipping_region,2);
[v,e,c] = polymesh_from_polysoup(p);
polymesh_save('voronoi_disk_40000',v,e,c);
