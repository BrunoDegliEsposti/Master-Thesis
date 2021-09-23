addpath('../Euler Equation');
addpath('../FVM Solver');
addpath('../Geometry');
addpath('../Meshes');

% Definizione dell'equazione differenziale
flux = @euler_flux2D;
t0 = 0;
freestream_u = u_from_rhovp([1.225,0,0,101325]);
mws = max_wave_speed2D(freestream_u);
T = 60/mws;
x0 = -20;
y0 = 0;
A = 0.1;
s = 0.2;
u0 = @(x,y) initial_conditions_u(x-x0,y-y0,freestream_u,A,s);

% Definizione del dominio discreto e dell'IVBP
[vertices,edges,cells] = polymesh_load('voronoi_ellipse_2500.mat');
cells.nu = 4;
cells.u = cell_integral(u0,cells.nu,vertices,edges,cells)./cells.area;
bc = {};
bc{1} = 'wall';

% Scelta dei metodi numerici
method.nq = 2;
method.order = 3;
method.reconstruction_strategy = @reconstruction_T1WENO3;
method.bc = bc;
method.flux = flux;
method.numerical_flux = @numerical_flux_rusanov;
method.ODE_solver = @SSPRK33;
method.courant_number = 1;

% Calcolo della soluzione numerica
prefix = 'ellipse-reflection';
tsnapshots = linspace(t0,T,121);
[vertices,edges,cells,niter] = solver(...
    t0,T,prefix,tsnapshots,vertices,edges,cells,method);



