addpath('../Geometry');

%% Griglia regolare
p = polysoup_from_grid(32,18,0,0,16,9);
polysoup_plot(p);

%% Corona circolare
p = polysoup_from_grid(60,10,0,0.5,2*pi,0.5);
p = polysoup_transform(p,@pol2cart);
polysoup_plot(p);

%% Disco
p = polysoup_from_grid(60,10,0,0,2*pi,1);
p = polysoup_transform(p,@pol2cart);
polysoup_plot(p);

%% Triangolazione arbitraria
TR = stlread('annulus.stl');
p = polysoup_from_triangulation(TR);
polysoup_plot(p);

%% Celle di Voronoi in un quadrato
N = 200;
hs = haltonset(2);
x = hs(1:N,1);
y = hs(1:N,2);
clipping_region = polyshape([0,1,1,0],[0,0,1,1]);
p = polysoup_from_voronoi(x,y,clipping_region);
polysoup_plot(p);
%hold on;
%scatter(x,y);

%% Regolarizzazione tramite algoritmo di Lloyd
N = 200;
hs = haltonset(2);
x = hs(1:N,1);
y = hs(1:N,2);
clipping_region = polyshape([0,1,1,0],[0,0,1,1]);
p = polysoup_from_voronoi_lloyd(x,y,clipping_region,3);
polysoup_plot(p);

%% Celle di Voronoi in un dominio arbitrario
N = 200;
hs = haltonset(2);
x = 2*hs(1:N,1)-1;
y = 2*hs(1:N,2)-1;
TR = stlread('blob.stl');
ps = polysoup_from_triangulation(TR);
clipping_region = polyshape_from_polysoup(ps);
p = polysoup_from_voronoi(x,y,clipping_region);
polysoup_plot(p);

%% Celle di Voronoi in un dominio arbitrario + LLoyd
N = 1000;
hs = haltonset(2);
x = 2*hs(1:N,1)-1;
y = 2*hs(1:N,2)-1;
TR = stlread('blob.stl');
ps = polysoup_from_triangulation(TR);
clipping_region = polyshape_from_polysoup(ps);
p = polysoup_from_voronoi_lloyd(x,y,clipping_region,4);
polysoup_plot(p);


