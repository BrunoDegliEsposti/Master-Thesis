function [q] = cell_area(vertices,edges,cells,i)
%CELL_AREA Calcola l'area delle celle con indici i
    if nargin < 4
        i = 1:cells.nc;
    end
    g = @(x,y) 1;
    q = cell_integral(g,1,vertices,edges,cells,i);
end
