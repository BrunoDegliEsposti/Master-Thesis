% Test a pagina 9 di Mengaldo et al.

addpath('../Euler Equation');
addpath('../FVM Solver');
addpath('../Geometry');
addpath('../Meshes');

% Definizione dell'equazione differenziale
flux = @euler_flux2D;
t0 = 0;
freestream_u = u_from_rhovMp([1.225,1,0,0.5,101325]);
mws = max_wave_speed2D(freestream_u);
T = 3/mws;
A = 0.1;
s = 0.02;
u0 = @(x,y) initial_conditions_u(x,y,freestream_u,A,s);

% Definizione del dominio discreto e dell'IVBP
[vertices,edges,cells] = polymesh_load('voronoi_disk_2500.mat');
% polymesh_plot(vertices,edges,cells);
cells.nu = 4;
cells.u = cell_integral(u0,cells.nu,vertices,edges,cells)./cells.area;
bc = containers.Map('KeyType','uint32','ValueType','any');
bc(1) = freestream_u;

% Scelta dei metodi numerici
edges.nq = 1;
edges = initialize_edge_quadrature(edges);
method.reconstruction_strategy = @reconstruction_constant;
method.bc = bc;
method.flux = flux;
method.numerical_flux = @numerical_flux_rusanov;
method.ODE_solver = @SSPRK11;
method.courant_number = 1;

% Calcolo della soluzione numerica
prefix = 'mengaldo1';
tsnapshots = linspace(t0,T,61);
[vertices,edges,cells,niter] = solver(...
    t0,T,prefix,tsnapshots,vertices,edges,cells,method);



