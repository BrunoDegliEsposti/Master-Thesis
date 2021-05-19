function [c] = speed_of_sound2D(u)
    %SPEED_OF_SOUND2D Calcola la velocità del suono in un gas nello stato u.
    p = pressure2D(u);
    c = realsqrt(adiabatic_index*p./u(:,1));
end
