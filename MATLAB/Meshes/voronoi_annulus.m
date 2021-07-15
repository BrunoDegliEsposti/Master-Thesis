%% Voronoi annulus
addpath('../Geometry');
cx = 0;
cy = 0;
ri = 0.1;
ro = 1;

N = 2500;
hs = haltonset(2);
x = (cx-ro) + 2*ro*hs(1:N,1);
y = (cy-ro) + 2*ro*hs(1:N,2);
pso = nsidedpoly(round(sqrt(N)*pi),'Center',[cx,cy],'Radius',ro);
psi = nsidedpoly(round(sqrt(N)*4*pi*ri/ro),'Center',[cx,cy],'Radius',ri);
clipping_region = subtract(pso,psi);
p = polysoup_from_voronoi_variable(x,y,clipping_region,12);
[v,e,c] = polymesh_from_polysoup(p);
polymesh_save('voronoi_annulus_2500',v,e,c);

N = 10000;
hs = haltonset(2);
x = (cx-ro) + 2*ro*hs(1:N,1);
y = (cy-ro) + 2*ro*hs(1:N,2);
pso = nsidedpoly(round(sqrt(N)*pi),'Center',[cx,cy],'Radius',ro);
psi = nsidedpoly(round(sqrt(N)*4*pi*ri/ro),'Center',[cx,cy],'Radius',ri);
clipping_region = subtract(pso,psi);
p = polysoup_from_voronoi_variable(x,y,clipping_region,12);
[v,e,c] = polymesh_from_polysoup(p);
polymesh_save('voronoi_annulus_10000',v,e,c);

N = 40000;
hs = haltonset(2);
x = (cx-ro) + 2*ro*hs(1:N,1);
y = (cy-ro) + 2*ro*hs(1:N,2);
pso = nsidedpoly(round(sqrt(N)*pi),'Center',[cx,cy],'Radius',ro);
psi = nsidedpoly(round(sqrt(N)*4*pi*ri/ro),'Center',[cx,cy],'Radius',ri);
clipping_region = subtract(pso,psi);
p = polysoup_from_voronoi_variable(x,y,clipping_region,12);
[v,e,c] = polymesh_from_polysoup(p);
polymesh_save('voronoi_annulus_40000',v,e,c);


