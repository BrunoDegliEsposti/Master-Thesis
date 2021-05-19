function [b] = ispolyconvex(ax,ay)
%ISPOLYCONVEX Controlla se un poligono in ingresso con vertici in ordine
% antiorario è convesso. Il formato è lo stesso di ispolycw().
    bx = [ax(2:end);ax(1)];
    by = [ay(2:end);ay(1)];
    cx = [ax(3:end);ax(1:2)];
    cy = [ay(3:end);ay(1:2)];
    b = all(counterclockwise(ax,ay,bx,by,cx,cy));
end
