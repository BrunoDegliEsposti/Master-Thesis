function cells = reconstruction_LLS3_initialize(vertices,edges,cells)
%RECONSTRUCTION_LLS3_INITIALIZE Inizializza il campo cells.camb nel caso
% di ricostruzioni di ordine 3
    cells.camb = zeros(cells.nc,3);
    gxx = @(x,y) x.*x;
    cells.camb(:,1) = cell_integral_mean(gxx,1,vertices,edges,cells);
    gxy = @(x,y) x.*y;
    cells.camb(:,2) = cell_integral_mean(gxy,1,vertices,edges,cells);
    gyy = @(x,y) y.*y;
    cells.camb(:,3) = cell_integral_mean(gyy,1,vertices,edges,cells);
end
