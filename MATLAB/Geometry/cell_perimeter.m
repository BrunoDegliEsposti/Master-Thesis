function [p] = cell_perimeter(vertices,edges,cells,i)
%CELL_PERIMETER Calcola il perimetro delle celle con indici i
    if nargin < 4
        i = 1:cells.nc;
    end
    p = zeros(length(i),1);
    for j=1:cells.mne
        e = abs(cells.e(i,j));
        mask = (e ~= 0);
        e = e(mask);
        l = edge_length(vertices,edges,e);
        p(mask) = p(mask) + l;
    end
end

