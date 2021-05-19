function [j] = edge_select_on_boundary(vertices,edges,f)
%EDGE_SELECT_ON_BOUNDARY Calcola tutti gli indici degli spigoli di bordo i cui
% estremi soddisfano entrambi la funzione booleana f: R^2 -> {true,false}.
    j = 1:edges.ne;
    [x1,y1] = edge_lerp(0,vertices,edges,j);
    [x2,y2] = edge_lerp(1,vertices,edges,j);
    mask1 = f(x1,y1);
    mask2 = f(x2,y2);
    j = j(mask1&mask2);
end
