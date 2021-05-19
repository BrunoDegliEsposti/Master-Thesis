function [q1,q2] = polymesh_quality(vertices,edges,cells)
%POLYMESH_QUALITY Calcola dei valori associati alla qualità della mesh
    q1 = min(2*cells.area./cells.perimeter);
    q2 = min((4*pi*cells.area)./(cells.perimeter.^2));
end
