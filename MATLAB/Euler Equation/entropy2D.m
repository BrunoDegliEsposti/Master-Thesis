function [s] = entropy2D(u)
%ENTROPY2D Calcola l'entropia specifica dell'aria nello stato u.
% Per convezione, l'aria nello stato p = 1 e rho = 1 ha entropia nulla.
    molar_mass_air = 2.896e-2; % kg mol^-1
    R = 8.3145; % J mol^-1 K^-1
    p = pressure2D(u); % Pa = kg m^-1 s^-2
    rho = u(:,1);
    s = (1/molar_mass_air) * (5/2) * R * reallog(p.*rho.^(-adiabatic_index));
end
