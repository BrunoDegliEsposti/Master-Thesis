function [rhovp] = rhovp_from_rhovMp(rhovMp)
%RHOVP_FROM_RHOVMP Calcola la velocità di un gas a partire dal numero di Mach.
    rho = rhovMp(:,1);
    vx = rhovMp(:,2);
    vy = rhovMp(:,3);
    M = rhovMp(:,4);
    p = rhovMp(:,5);
    rhovp = zeros(size(rho,1),4);
    c = realsqrt(adiabatic_index*p./rho);
    vnorm = hypot(vx,vy);
    rhovp(:,1) = rho;
    rhovp(:,2) = (vx./vnorm).*c.*M;
    rhovp(:,3) = (vy./vnorm).*c.*M;
    rhovp(:,4) = p;
end
