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
[vertices,edges,cells] = polymesh_load('regular_square_200x200.mat');
cells.nu = 4;
cells.u = cell_integral(u0,cells.nu,vertices,edges,cells)./cells.area;
bc = containers.Map('KeyType','uint32','ValueType','any');
bc(1) = u_exact;

% Scelta dei metodi numerici
edges.nq = 1;
edges = initialize_edge_quadrature(edges);
L = @constant_WR_Rusanov_FVM;
%cells = initialise_cell_stencils(vertices,edges,cells); TODO...
ODE_solver = @SSPRK11;
courant_number = 1;

% Calcolo della soluzione numerica
prefix = 'vortex-grid';
tsnapshots = linspace(t0,T,11);
[vertices,edges,cells,niter] = solver(t0,T,prefix,tsnapshots,...
    vertices,edges,cells,ODE_solver,courant_number,L,bc,flux);

% Stima dell'errore L1
g = @(x,y,ubar) abs(u_exact(x,y,T)-ubar);
errL1 = sum(cell_integral_u(g,cells.nu,vertices,edges,cells));

% Stima dell'errore L2
g = @(x,y,ubar) (u_exact(x,y,T)-ubar).^2;
errL2 = realsqrt(sum(cell_integral_u(g,cells.nu,vertices,edges,cells)));



