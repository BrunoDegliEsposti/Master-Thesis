function [rhovp] = rhovp_from_u(u)
%RHOVP_FROM_U Calcola le variabili primitive (densità, velocità, pressione)
% a partire dalle variabili conservate (densità, quantità di moto, energia totale).
    rhovp = zeros(size(u));
    rhovp(:,1) = u(:,1);
    rhovp(:,2) = u(:,2)./u(:,1);
    rhovp(:,3) = u(:,3)./u(:,1);
    rhovp(:,4) = pressure2D(u);
end
