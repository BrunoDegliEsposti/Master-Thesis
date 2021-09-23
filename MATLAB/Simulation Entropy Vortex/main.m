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
cx0 = 0;
cy0 = 0;
cvx = 1;
cvy = 0;
beta = 5;
u_exact = @(x,y,t) vortex(x,y,t,cx0,cy0,cvx,cvy,beta);
u0 = @(x,y) u_exact(x,y,t0);

% Definizione del dominio discreto e dell'IVBP
[vertices,edges,cells] = polymesh_load('regular_square_100x100.mat');
cells.nu = 4;
cells.u = cell_integral_mean(u0,cells.nu,vertices,edges,cells);
bc = {};
bc{1} = u_exact;

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
prefix = 'vortex-grid';
tsnapshots = linspace(t0,T,11);
[vertices,edges,cells,niter] = solver(...
    t0,T,prefix,tsnapshots,vertices,edges,cells,method);

% Stima degli errori L1, L2 e Linf
uT = @(x,y) u_exact(x,y,T);
errL1 = cell_norm_L1(...
    cells.u - cell_integral_mean(uT,cells.nu,vertices,edges,cells), cells);
errL2 = cell_norm_L2(...
    cells.u - cell_integral_mean(uT,cells.nu,vertices,edges,cells), cells);
errLinf = cell_norm_Linf(...
    cells.u - cell_integral_mean(uT,cells.nu,vertices,edges,cells));




