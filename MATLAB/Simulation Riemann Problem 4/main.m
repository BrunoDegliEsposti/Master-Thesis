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
polysoup = polysoup_from_grid(200,200,-1,-1,3,3);
[vertices,edges,cells] = polymesh_from_polysoup(polysoup);
cells.nu = 4;
cells.u = cell_integral(u0,cells.nu,vertices,edges,cells)./cells.area;
bc = {};
bc{1} = 'absorbing';

% Scelta dei metodi numerici
edges.nq = 1;
edges = initialize_edge_quadrature(edges);
%cells = reconstruction_LLS3_initialize(vertices,edges,cells);
method.reconstruction_strategy = @reconstruction_LLS2P;
method.bc = bc;
method.flux = flux;
method.numerical_flux = @numerical_flux_rusanov;
method.ODE_solver = @SSPRK22;
method.courant_number = 1;

% Calcolo della soluzione numerica
prefix = 'riemann-problem-4';
tsnapshots = linspace(t0,T,101);
[vertices,edges,cells,niter] = solver(...
    t0,T,prefix,tsnapshots,vertices,edges,cells,method);




