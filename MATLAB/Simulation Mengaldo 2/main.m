% Test a pagina 12 di Mengaldo et al.

addpath('../Euler Equation');
addpath('../FVM Solver');
addpath('../Geometry');
addpath('../Meshes');

% Definizione dell'equazione differenziale
flux = @euler_flux2D;
t0 = 0;
freestream_u = u_from_rhovMp([1.225,1,0,0.1,101325]);
mws = max_wave_speed2D(freestream_u);
T = 5/mws;
u0 = @(x,y) freestream_u;

% Definizione del dominio discreto e dell'IVBP
[vertices,edges,cells] = polymesh_load('voronoi_annulus_2500.mat');
j = edge_select_on_boundary(vertices,edges, @(x,y) hypot(x,y)<0.2 );
edges.type(j) = 2;
% polymesh_plot(vertices,edges,cells);
cells.nu = 4;
cells.u = cell_integral(u0,cells.nu,vertices,edges,cells)./cells.area;
bc = containers.Map('KeyType','uint32','ValueType','any');
bc(1) = freestream_u;
bc(2) = 'wall';

% Scelta dei metodi numerici
edges.nq = 1;
edges = initialize_edge_quadrature(edges);
L = @constant_WR_Rusanov_FVM;
%cells = initialise_cell_stencils(vertices,edges,cells); TODO...
ODE_solver = @SSPRK11;
courant_number = 1;

% Calcolo della soluzione numerica
prefix = 'mengaldo2';
tsnapshots = linspace(t0,T,11);
[vertices,edges,cells,niter] = solver(t0,T,prefix,tsnapshots,...
    vertices,edges,cells,ODE_solver,courant_number,L,bc,flux);



