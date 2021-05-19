function [nx,ny] = edge_normal(vertices,edges,j)
%EDGE_NORMAL Calcola la normale degli spigoli con indici j.
    if nargin < 3
        j = 1:edges.ne;
    end
    v1 = edges.v1(j);
    v2 = edges.v2(j);
    v1x = vertices.x(v1);
    v1y = vertices.y(v1);
    v2x = vertices.x(v2);
    v2y = vertices.y(v2);
    nx = -(v2y-v1y);
    ny =   v2x-v1x;
    nl = realsqrt(nx.*nx+ny.*ny);
    nx = nx./nl;
    ny = ny./nl;
end
