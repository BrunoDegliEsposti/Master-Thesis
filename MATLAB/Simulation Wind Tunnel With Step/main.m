addpath('../Euler Equation');
addpath('../FVM Solver');
addpath('../Geometry');
addpath('../Meshes');

% Definizione dell'equazione differenziale
flux = @euler_flux2D;
t0 = 0;
freestream_u = u_from_rhovp([adiabatic_index,3,0,1]);
T = 4;
u0 = @(x,y) freestream_u;

% Definizione del dominio discreto e dell'IVBP
[vertices,edges,cells] = polymesh_load('voronoi_wind_tunnel_with_step_10000.mat');
j = edge_select_on_boundary(vertices,edges, @(x,y) x < 1e-8);
edges.type(j) = 2;
j = edge_select_on_boundary(vertices,edges, @(x,y) x > 3-1e-8);
edges.type(j) = 3;
%polymesh_plot(vertices,edges,cells);
cells.nu = 4;
cells.u = cell_integral(u0,cells.nu,vertices,edges,cells)./cells.area;
bc = {};
bc{1} = 'wall';
bc{2} = freestream_u;
bc{3} = 'absorbing';

% Scelta dei metodi numerici
method.nq = 1;
method.order = 2;
method.reconstruction_strategy = @reconstruction_T1WENO;
method.WENO_epsilon = 1e-6*max(cells.h)^2;
method.WENO_power = 4;
method.bc = bc;
method.flux = flux;
method.numerical_flux = @numerical_flux_rusanov;
method.ODE_solver = @SSPRK22;
method.courant_number = 1;

% Calcolo della soluzione numerica
prefix = 'wind-tunnel-with-step';
tsnapshots = linspace(t0,T,41);
[vertices,edges,cells,niter] = solver(...
    t0,T,prefix,tsnapshots,vertices,edges,cells,method);



