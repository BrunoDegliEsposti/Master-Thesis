function [b] = edge_plusside(cx,cy,vertices,edges,j)
%EDGE_PLUSSIDE Determina se i punti (cx,cy) stanno dal lato uplus
% degli spigoli con indici j. L'esito è un booleano.
    if nargin < 5
        j = 1:edges.ne;
    end
    assert(length(cx)==length(j));
    assert(length(cy)==length(j));
    v1 = edges.v1(j);
    v2 = edges.v2(j);
    v1x = vertices.x(v1);
    v1y = vertices.y(v1);
    v2x = vertices.x(v2);
    v2y = vertices.y(v2);
    b = counterclockwise(cx,cy,v1x,v1y,v2x,v2y);
end
