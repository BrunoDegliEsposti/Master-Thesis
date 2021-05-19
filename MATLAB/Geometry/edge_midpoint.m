function [mx,my] = edge_midpoint(vertices,edges,j)
%EDGE_MIDPOINT Calcola il punto medio degli spigoli con indici j.
    if nargin < 3
        j = 1:edges.ne;
    end
    [mx,my] = edge_lerp(0.5,vertices,edges,j);
end
