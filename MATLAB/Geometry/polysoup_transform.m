function [p] = polysoup_transform(p,f,tol)
%POLYSOUP_TRANSFORM Trasforma i vertici di una polygonal soup tramite una
% funzione f(x,y) -> \R^2 che accetta vettori colonna e restituisce una
% comma-separated list di vettori colonna.
% I vertici che dopo la trasformazione si trovano a meno di tol
% l'uno dall'altro vengono fusi in un unico vertice.
    if nargin < 3
        tol = 1e-8;
    end
    [p.vx,p.vy] = f(p.vx,p.vy);
    [p.cx,p.cy] = f(p.cx,p.cy);
    p = polysoup_merge_vertices(p,tol);
    p = polysoup_fix_CCW(p);
end
