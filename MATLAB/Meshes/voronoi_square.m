%% Voronoi square
addpath('../Geometry');
blcx = -5;
blcy = -5;
Lx = 10;
Ly = 10;
clipping_region = polyshape(...
    [blcx,blcx+Lx,blcx+Lx,blcx],[blcy,blcy,blcy+Ly,blcy+Ly]);
hs = haltonset(2);

N = 2500;
x = blcx + Lx*hs(1:N,1);
y = blcy + Ly*hs(1:N,2);
p = polysoup_from_voronoi_uniform(x,y,clipping_region,2);
[v,e,c] = polymesh_from_polysoup(p);
polymesh_save('voronoi_square_2500.mat',v,e,c);

N = 10000;
x = blcx + Lx*hs(1:N,1);
y = blcy + Ly*hs(1:N,2);
p = polysoup_from_voronoi_uniform(x,y,clipping_region,2);
[v,e,c] = polymesh_from_polysoup(p);
polymesh_save('voronoi_square_10000.mat',v,e,c);

N = 40000;
x = blcx + Lx*hs(1:N,1);
y = blcy + Ly*hs(1:N,2);
p = polysoup_from_voronoi_uniform(x,y,clipping_region,2);
[v,e,c] = polymesh_from_polysoup(p);
polymesh_save('voronoi_square_40000.mat',v,e,c);

N = 160000;
x = blcx + Lx*hs(1:N,1);
y = blcy + Ly*hs(1:N,2);
p = polysoup_from_voronoi_uniform(x,y,clipping_region,2);
[v,e,c] = polymesh_from_polysoup(p);
polymesh_save('voronoi_square_160000.mat',v,e,c);

