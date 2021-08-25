% Soluzione dell'equazione di Eulero per gas perfetti in 2D
% su dominio rettangolare discretizzato con una griglia regolare uniforme
% e con condizioni iniziali/al bordo corrispondenti a vortice isoentropico.
% La parte spaziale dell'equazione viene discretizzata con il metodo
% dei volumi finiti basato su ricostruzioni ENO/WENO e flusso di Rusanov.
% La parte temporale dell'equazione viene discretizzata con
% un metodo di Runge-Kutta. L'intervallo di simulazione è [0,T].

addpath('../Euler Equation');
addpath('../FVM Solver');
addpath('../Geometry');
addpath('../Meshes');

% Definizione dell'equazione differenziale
flux = @euler_flux2D;
t0 = 0;
T = 1;
rho_freestream = 1;
vx_freestream = 1;
vy_freestream = 1;
p_freestream = 1;
amplitude = 0.01;
u_exact = @(x,y,t) entropy_waves(x,y,t,rho_freestream,...
    vx_freestream,vy_freestream,p_freestream,amplitude);
u0 = @(x,y) entropy_waves(x,y,t0,rho_freestream,...
    vx_freestream,vy_freestream,p_freestream,amplitude);

% Definizione del dominio discreto e dell'IVBP
[vertices,edges,cells] = polymesh_load('regular_square_50x50.mat');
cells.nu = 4;
cells.u = cell_integral_mean(u0,cells.nu,vertices,edges,cells);
bc = {};
bc{1} = u_exact;

% Scelta dei metodi numerici
edges.nq = 1;
edges = initialize_edge_quadrature(edges);
method.reconstruction_strategy = @reconstruction_LLS2;
method.bc = bc;
method.flux = flux;
method.numerical_flux = @numerical_flux_rusanov;
method.ODE_solver = @SSPRK22;
method.courant_number = 1;

% Calcolo della soluzione numerica
prefix = 'vortex-grid';
tsnapshots = linspace(t0,T,101);
[vertices,edges,cells,niter] = solver(...
    t0,T,prefix,tsnapshots,vertices,edges,cells,method);

% Stima degli errori L1 e Linf
uT = @(x,y) u_exact(x,y,T);
errL1 = cell_norm_L1(...
    cells.u - cell_integral_mean(uT,cells.nu,vertices,edges,cells), cells);
errL2 = cell_norm_L2(...
    cells.u - cell_integral_mean(uT,cells.nu,vertices,edges,cells), cells);
errLinf = cell_norm_Linf(...
    cells.u - cell_integral_mean(uT,cells.nu,vertices,edges,cells));



