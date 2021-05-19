function [l] = edge_length(vertices,edges,j)
%EDGE_LENGTH Calcola la lunghezza degli spigoli con indici j.
    if nargin < 3
        j = 1:edges.ne;
    end
    v1 = edges.v1(j);
    v2 = edges.v2(j);
    v1x = vertices.x(v1);
    v1y = vertices.y(v1);
    v2x = vertices.x(v2);
    v2y = vertices.y(v2);
    l = realsqrt((v1x-v2x).^2+(v1y-v2y).^2);
end
