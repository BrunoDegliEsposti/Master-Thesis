function [gbar] = cell_integral_mean(g,ndim,vertices,edges,cells,i)
%CELL_INTEGRAL_MEAN Calcola la media integrale di una funzione vettoriale
% g(x,y) di ndim componenti sulle celle con indici i (dove "i" può essere
% un vettore). Le celle devono essere poligoni convessi e devono contenere
% sul bordo o al loro interno il punto (cells.cx,cells.cy).
% Per i dettagli, vedi cell_integral.m
    if nargin < 6
        i = 1:cells.nc;
    end
    gbar = cell_integral(g,ndim,vertices,edges,cells,i) ./ cells.area(i);
end
