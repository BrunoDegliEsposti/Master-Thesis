addpath('../Euler Equation');
addpath('../FVM Solver');
addpath('../Geometry');
addpath('../Meshes');

% Errore locale di troncamento al tempo t0
flux = @euler_flux2D;
t0 = 0;
cx0 = 0;
cy0 = 0;
cvx = 1;
cvy = 0;
beta = 5;
u_exact = @(x,y,t) vortex(x,y,t,cx0,cy0,cvx,cvy,beta);
u0 = @(x,y) u_exact(x,y,t0);

% Definizione della mesh, dei punti di quadratura e delle medie u
[vertices,edges,cells] = polymesh_load('regular_square_400x400.mat');
cells.nu = 4;
cells.u = cell_integral_mean(u0,cells.nu,vertices,edges,cells);
edges.nq = 3;
edges = initialize_edge_quadrature(edges);

% Calcolo delle ricostruzioni interne e sul bordo
cells = interpolation_linear_initialize(vertices,edges,cells);
[vertices,edges,cells] = interpolation_linear(vertices,edges,cells);
for j = edges.nie+(1:edges.nbe)
    for k = 1:edges.nq
        [x,y] = edge_lerp(edges.qx(k),vertices,edges,j);
        edges.um(j,:,k) = u0(x,y);
    end
end

% Approssimazione dell'integrale del flusso su ogni interfaccia tra celle
edges.tnf = zeros(edges.ne,cells.nu);
for k = 1:edges.nq
    [nf,~] = numerical_flux_rusanov(...
        flux,edges.um(:,:,k),edges.up(:,:,k),edges.nx,edges.ny);
    edges.tnf = edges.tnf + edges.qw(k) * nf;
end
edges.tnf = edges.length .* edges.tnf;

% Calcolo del massimo errore locale di troncamento al variare delle celle
taumax = zeros(1,cells.nu);
imax = 0;
for i = 1:cells.nc
    tau = zeros(1,cells.nu);
    for j = 1:cells.ne(i)
        e = cells.e(i,j);
        if e > 0
            % normale entrante nella cella i (l'opposto del teorema della div.)
            tau = tau - edges.tnf(e);
        elseif e < 0
            % normale uscente dalla cella i (come nel teorema della div.)
            tau = tau + edges.tnf(-e);
        end
        for k = 1:edges.nq
            if e > 0
                % normale entrante nella cella i (l'opposto del teorema della div.)
                [x,y] = edge_lerp(edges.qx(k),vertices,edges,e);
                F = flux(u0(x,y));
                Fnormal = F(1,:,1)*(-edges.nx(e)) + F(1,:,2)*(-edges.ny(e));
                tau = tau - edges.length(e) * edges.qw(k) * Fnormal;
            elseif e < 0
                % normale uscente dalla cella i (come nel teorema della div.)
                [x,y] = edge_lerp(edges.qx(k),vertices,edges,-e);
                F = flux(u0(x,y));
                Fnormal = F(1,:,1)*edges.nx(-e) + F(1,:,2)*edges.ny(-e);
                tau = tau - edges.length(-e) * edges.qw(k) * Fnormal;
            end
        end
    end
    tau = (1/cells.area(i)) * abs(tau);
    if norm(tau) > norm(taumax)
        taumax = tau;
        imax = i;
    end
end







