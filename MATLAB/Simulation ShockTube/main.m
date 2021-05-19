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
T = 0.2;
u0 = @(x,y) [1,   0, 0,   1/(adiabatic_index-1)].*(x<0.5) + ...
            [1/8, 0, 0, 0.1/(adiabatic_index-1)].*(x>=0.5);

% Definizione del dominio discreto e dell'IVBP
polysoup = polysoup_from_grid(1024,3,0,0,1,3/1024);
[vertices,edges,cells] = polymesh_from_polysoup(polysoup);
cells.nu = 4;
cells.u = cell_integral(u0,cells.nu,vertices,edges,cells)./cells.area;
bc = containers.Map('KeyType','uint32','ValueType','any');
bc(1) = 'absorbing';

% Scelta dei metodi numerici
edges.nq = 1;
edges = initialize_edge_quadrature(edges);
L = @constant_WR_Rusanov_FVM;
%cells = initialise_cell_stencils(vertices,edges,cells); TODO...
ODE_solver = @SSPRK11;
courant_number = 1;

% Calcolo della soluzione numerica
prefix = 'shock-tube';
tsnapshots = linspace(t0,T,11);
[vertices,edges,cells,niter] = solver(t0,T,prefix,tsnapshots,...
    vertices,edges,cells,ODE_solver,courant_number,L,bc,flux);




