function [u] = u_from_TvMp(TvMp)
%U_FROM_RHOVMP Calcola le variabili conservate (densità, quantità di moto,
% energia totale) a partire da temperatura (T), direzione della velocità (v),
% numero di Mach (M) e pressione (p).
    T = TvMp(:,1);
    vx = TvMp(:,2);
    vy = TvMp(:,3);
    M = TvMp(:,4);
    p = TvMp(:,5);
    u = zeros(size(T,1),4);
    molar_mass_air = 2.896e-2; % kg mol^-1
    R = 8.3145; % J mol^-1 K^-1
    rho = molar_mass_air * (1/R) * p ./ T;
    c = realsqrt(adiabatic_index*p./rho);
    vnorm = hypot(vx,vy);
    vx = (vx./vnorm).*c.*M;
    vy = (vy./vnorm).*c.*M;
    u(:,1) = rho;
    u(:,2) = rho.*vx;
    u(:,3) = rho.*vy;
    u(:,4) = p/(adiabatic_index-1) + 0.5*rho.*(vx.*vx+vy.*vy);
end

