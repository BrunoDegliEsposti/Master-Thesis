function [u] = u_from_rhovMp(rhovMp)
%U_FROM_RHOVMP Calcola le variabili conservate (densità, quantità di moto,
% energia totale) a partire da densità (rho), direzione della velocità (v),
% numero di Mach (M) e pressione (p).
    rho = rhovMp(:,1);
    vx = rhovMp(:,2);
    vy = rhovMp(:,3);
    M = rhovMp(:,4);
    p = rhovMp(:,5);
    u = zeros(size(rho,1),4);
    c = realsqrt(adiabatic_index*p./rho);
    vnorm = hypot(vx,vy);
    vx = (vx./vnorm).*c.*M;
    vy = (vy./vnorm).*c.*M;
    u(:,1) = rho;
    u(:,2) = rho.*vx;
    u(:,3) = rho.*vy;
    u(:,4) = p/(adiabatic_index-1) + 0.5*rho.*(vx.*vx+vy.*vy);
end

