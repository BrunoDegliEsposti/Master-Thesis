function [M] = mach_number2D(u)
%MACH_NUMBER2D Calcola il numero di Mach di un gas nello stato u.
    c = speed_of_sound2D(u);
    vx = u(:,2)./u(:,1);
    vy = u(:,3)./u(:,1);
    v = realsqrt(vx.^2+vy.^2);
    M = v./c;
end

