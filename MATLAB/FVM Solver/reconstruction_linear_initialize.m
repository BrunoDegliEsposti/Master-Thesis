function cells = reconstruction_linear_initialize(vertices,edges,cells)
%RECONSTRUCTION_LINEAR_INITIALIZE Inizializza il campo cells.imm nel caso
% di ricostruzioni lineari (grado 1, quindi ordine 2).
    cells.imm = zeros(cells.nc,3);
    cells.imm(:,1) = 1;
    gx = @(x,y) x;
    cells.imm(:,2) = cell_integral_mean(gx,1,vertices,edges,cells);
    gy = @(x,y) y;
    cells.imm(:,3) = cell_integral_mean(gy,1,vertices,edges,cells);
end
