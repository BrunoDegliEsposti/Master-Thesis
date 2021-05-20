function [err] = cell_norm_Linf(u)
%CELL_NORM_LINF Calcola la norma L^\infty discreta dei valori u (che vanno
% pensati come medie integrali). Il numero di colonne di err è uguale al
% numero di colonne di u.
    err = max(abs(u));
end
