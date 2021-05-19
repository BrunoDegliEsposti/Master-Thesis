function [mws] = max_wave_speed2D(u)
    %MAX_WAVE_SPEED2D Calcola la massima velocità a cui si propaga
    % un gas nello stato u.
    vx = u(:,2)./u(:,1);
    vy = u(:,3)./u(:,1);
    c = speed_of_sound2D(u);
    mws = realsqrt(vx.*vx+vy.*vy)+c;
end
