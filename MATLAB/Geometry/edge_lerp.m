function [x,y] = edge_lerp(t,vertices,edges,j)
%EDGE_LERP Calcola dei punti sugli spigoli con indici j.
% Ogni punto è la combinazione convessa con parametro t dei due vertici.
% t = 0 corrisponde a v1, t = 1 corrisponde a v2.
    if nargin < 4
        j = 1:edges.ne;
    end
    v1 = edges.v1(j);
    v2 = edges.v2(j);
    v1x = vertices.x(v1);
    v1y = vertices.y(v1);
    v2x = vertices.x(v2);
    v2y = vertices.y(v2);
    x = (1-t).*v1x + t.*v2x;
    y = (1-t).*v1y + t.*v2y;
end
