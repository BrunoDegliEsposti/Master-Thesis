%% Voronoi ellipse
% terna pitagorica 20 (fuochi), 21 (altezza), 29 (larghezza)
addpath('../Geometry');
N = 10000;
hs = haltonset(2);
x = 60*hs(1:N,1)-30;
y = 60*hs(1:N,2)-30;
TR = stlread('ellipse.stl');
ps = polysoup_from_triangulation(TR);
clipping_region = polyshape_from_polysoup(ps);
p = polysoup_from_voronoi_uniform(x,y,clipping_region,4);
[v,e,c] = polymesh_from_polysoup(p);
polymesh_save('voronoi_ellipse_10000',v,e,c);
polymesh_plot(v,e,c);
