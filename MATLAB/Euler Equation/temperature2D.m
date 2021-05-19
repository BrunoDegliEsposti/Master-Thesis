function [T] = temperature2D(u)
    %TEMPERATURE2D Calcola la temperatura di un gas nello stato u.
    molar_mass_air = 2.896e-2; % kg mol^-1
    R = 8.3145; % J mol^-1 K^-1
    p = pressure2D(u); % Pa = kg m^-1 s^-2
    T = molar_mass_air * (1/R) * p ./ u(:,1);
end
