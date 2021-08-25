%% Voronoi wind tunnel with step
addpath('../Geometry');
N = 10000;
hs = haltonset(2);
x = 3*hs(1:N,1);
y = 3*hs(1:N,2)-1;
TR = stlread('wind_tunnel_with_step.stl');
ps = polysoup_from_triangulation(TR);
clipping_region = polyshape_from_polysoup(ps);
p = polysoup_from_voronoi_uniform(x,y,clipping_region,4);
[v,e,c] = polymesh_from_polysoup(p);
polymesh_save('voronoi_wind_tunnel_with_step_10000',v,e,c);
polymesh_plot(v,e,c);
