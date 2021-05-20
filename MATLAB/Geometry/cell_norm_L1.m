function [err] = cell_norm_L1(u,cells,i)
%CELL_NORM_L1 Calcola la norma L1 discreta dei valori u (che vanno pensati
% come medie integrali) sulle celle con indici i (dove "i" può essere
% un vettore). Il numero di colonne di err è uguale al numero di colonne di u.
    if nargin < 3
        i = 1:cells.nc;
    end
    assert(length(i) == size(u,1));
    err = sum(cells.area(i) .* abs(u)) / sum(cells.area(i));
end
