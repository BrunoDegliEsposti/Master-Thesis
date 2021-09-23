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
T = 0.5;
v_shock = 2;

g = adiabatic_index;
p_left = (1 + (v_shock^2-1)*(2*g)/(g+1))/g;
rho_left = (1 + p_left*(g+1)/(g-1)) / (p_left + (g+1)/(g-1));
v_left = (g*p_left - 1) / (g*realsqrt(1 + (g*p_left-1)*(g+1)/(2*g)));
u_left = u_from_rhovp([rho_left,v_left,0,p_left]);

u_right = u_from_rhovp([1,0,0,1/adiabatic_index]);
u0 = @(x,y) u_left.*(x<0) + u_right.*(x>=0);

% Definizione del dominio discreto e dell'IVBP
polysoup = polysoup_from_grid(128,64,-1,-0.5,2,1);
[vertices,edges,cells] = polymesh_from_polysoup(polysoup);
j = edge_select_on_boundary(vertices,edges, @(x,y) x < -1+1e-8);
edges.type(j) = 2;
j = edge_select_on_boundary(vertices,edges, @(x,y) x > 1-1e-8);
edges.type(j) = 3;
%polymesh_plot(vertices,edges,cells);
cells.nu = 4;
cells.u = cell_integral(u0,cells.nu,vertices,edges,cells)./cells.area;
bc = {};
bc{1} = 'wall';
bc{2} = u_left;
bc{3} = u_right;

% Scelta dei metodi numerici
method.nq = 1;
method.order = 1;
method.reconstruction_strategy = @reconstruction_LLS1;
method.bc = bc;
method.flux = flux;
method.numerical_flux = @numerical_flux_rusanov;
method.ODE_solver = @SSPRK11;
method.courant_number = 1;

% Calcolo della soluzione numerica
prefix = 'propagating-normal-shock';
tsnapshots = linspace(t0,T,51);
[vertices,edges,cells,niter] = solver(...
    t0,T,prefix,tsnapshots,vertices,edges,cells,method);




