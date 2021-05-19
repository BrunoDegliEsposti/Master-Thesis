function [u] = u_from_rhovp(rhovp)
%U_FROM_RHOVP Calcola le variabili conservate (densità, quantità di moto,
% energia totale) a partire da densità (rho), velocità (v) e pressione (p).
    rho = rhovp(:,1);
    vx = rhovp(:,2);
    vy = rhovp(:,3);
    p = rhovp(:,4);
    u = zeros(size(rhovp));
    u(:,1) = rho;
    u(:,2) = rho.*vx;
    u(:,3) = rho.*vy;
    u(:,4) = p/(adiabatic_index-1) + 0.5*rho.*(vx.*vx+vy.*vy);
end
