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
T = 0.25;
u0 = @(x,y) ...
    u_from_rhovp([1.1,    0.8939, 0.8939, 1.1 ]).*(x< 0.5).*(y< 0.5) + ...
    u_from_rhovp([0.5065, 0,      0.8939, 0.35]).*(x>=0.5).*(y< 0.5) + ...
    u_from_rhovp([0.5065, 0.8939, 0,      0.35]).*(x< 0.5).*(y>=0.5) + ...
    u_from_rhovp([1.1,    0,      0,      1.1 ]).*(x>=0.5).*(y>=0.5);

% Definizione del dominio discreto e dell'IVBP
polysoup = polysoup_from_grid(400,400,-0.5,-0.5,2,2);
[vertices,edges,cells] = polymesh_from_polysoup(polysoup);
cells.nu = 4;
cells.u = cell_integral(u0,cells.nu,vertices,edges,cells)./cells.area;
bc = {};
bc{1} = 'absorbing';

% Scelta dei metodi numerici
method.nq = 1;
method.order = 2;
method.reconstruction_strategy = @reconstruction_T1WENO;
%method.least_squares_type = 'P';
method.WENO_epsilon = 1e-6 * max(cells.h)^2;
method.WENO_power = 4;
method.bc = bc;
method.flux = flux;
method.numerical_flux = @numerical_flux_rusanov;
method.ODE_solver = @SSPRK22;
method.courant_number = 1;

% Calcolo della soluzione numerica
prefix = 'riemann-problem-4';
tsnapshots = linspace(t0,T,26);
[vertices,edges,cells,niter] = solver(...
    t0,T,prefix,tsnapshots,vertices,edges,cells,method);




