addpath('../Geometry');

%% Griglia regolare
p = polysoup_from_grid(32,18,0,0,16,9);
[v,e,c] = polymesh_from_polysoup(p);
j = edge_select_on_boundary(v,e, @(x,y) y>9-1e-8 );
e.type(j) = 2;
j = edge_select_on_boundary(v,e, @(x,y) x>16-1e-8 );
e.type(j) = 3;
j = edge_select_on_boundary(v,e, @(x,y) y<1e-8 );
e.type(j) = 4;
polymesh_save('grid.mat',v,e,c);
polymesh_plot(v,e,c);

%% Caricamento da file
[v,e,c] = polymesh_load('grid.mat');
polymesh_plot(v,e,c);

%% Corona circolare
p = polysoup_from_grid(60,10,0,0.5,2*pi,0.5);
p = polysoup_transform(p,@pol2cart);
[v,e,c] = polymesh_from_polysoup(p);
j = edge_select_on_boundary(v,e, @(x,y) hypot(x,y)<0.6 );
e.type(j) = 2;
polymesh_plot(v,e,c);

%% Disco
p = polysoup_from_grid(10,4,0,0,2*pi,1);
p = polysoup_transform(p,@pol2cart);
g = @(x,y) deal(x+1,y);
p = polysoup_transform(p,g);
[v,e,c] = polymesh_from_polysoup(p);
polymesh_plot(v,e,c);

%% Celle di Voronoi in un quadrato + LLoyd
N = 200;
hs = haltonset(2);
x = hs(1:N,1);
y = hs(1:N,2);
clipping_region = polyshape([0,1,1,0],[0,0,1,1]);
p = polysoup_from_voronoi_lloyd(x,y,clipping_region,3);
[v,e,c] = polymesh_from_polysoup(p);
polymesh_plot(v,e,c);




