function [u] = vortex(x,y,t,cx0,cy0,cvx,cvy,beta)
    %VORTEX Vortice isoentropico
    % Soluzione di riferimento per l'equazione di Eulero 2D.
    % Il centro è inizialmente in (cx0,cy0) e si muove a velocità (cvx,cvy).
    % Il coefficiente beta è associato all'intensità del vortice.
    % Valori di beta troppo grandi formano il vuoto, quindi vanno evitati.
    u = zeros(length(x),4);
    % Coordinate x e y nel SDR solidale al centro del vortice
    xc = x-cvx*t-cx0;
    yc = y-cvy*t-cy0;
    % Distanza dal centro del vortice
    r = realsqrt(xc.^2+yc.^2);
    % Pressione
    g = (adiabatic_index-1)/(16*adiabatic_index*pi^2);
    u(:,1) = (1 - g * beta^2 * exp(2*(1-r.^2))) .^ (1/(adiabatic_index-1));
    % Velocità e quantità di moto
    vx = cvx - beta*exp(1-r.^2).*yc/(2*pi);
    vy = cvy + beta*exp(1-r.^2).*xc/(2*pi);
    u(:,2) = u(:,1).*vx;
    u(:,3) = u(:,1).*vy;
    % Pressione ed energia totale specifica
    p = u(:,1).^adiabatic_index;
    u(:,4) = p/(adiabatic_index-1) + (u(:,2).*u(:,2)+u(:,3).*u(:,3))./(2*u(:,1));
    % Controllo di correttezza
    check_state2D(u);
end
