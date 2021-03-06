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
T = 10/mws;
u0 = @(x,y) freestream_u;

% Definizione del dominio discreto e dell'IVBP
[vertices,edges,cells] = polymesh_load('voronoi_annulus_2500.mat');
j = edge_select_on_boundary(vertices,edges, @(x,y) hypot(x,y)<0.2 );
edges.type(j) = 2;
% polymesh_plot(vertices,edges,cells);
cells.nu = 4;
cells.u = cell_integral(u0,cells.nu,vertices,edges,cells)./cells.area;
bc = {};
bc{1} = freestream_u;
bc{2} = 'wall';

% Scelta dei metodi numerici
method.nq = 1;
method.order = 2;
method.reconstruction_strategy = @reconstruction_LLS2;
method.least_squares_type = 'P';
method.bc = bc;
method.flux = flux;
method.numerical_flux = @numerical_flux_rusanov;
method.ODE_solver = @SSPRK22;
method.courant_number = 1;

% Calcolo della soluzione numerica
prefix = 'mengaldo2';
tsnapshots = linspace(t0,T,101);
[vertices,edges,cells,niter] = solver(...
    t0,T,prefix,tsnapshots,vertices,edges,cells,method);



