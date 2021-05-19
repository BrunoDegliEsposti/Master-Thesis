function [q_polygon] = cell_integral(g,ndim,vertices,edges,cells,i)
%CELL_INTEGRAL Calcola l'integrale di una funzione vettoriale g(x,y)
% di ndim componenti sulle celle con indici i (dove "i" può essere un vettore).
% Le celle devono essere poligoni convessi e devono contenere sul bordo o
% al loro interno il punto (cells.cx,cells.cy).
% Il risultato q_polygon è il valore dell'integrale su ogni cella;
% questi integrali locali possono poi essere sommati per ottenere
% un integrale globale.
% L'integrale su ogni cella poligonale viene calcolato tramite
% una triangolazione "a ventaglio" del poligono intorno al centro,
% seguita da formule di quadratura su ogni triangolo del ventaglio.
% La formula per triangoli, dovuta a Radon, usa 7 punti ed è esatta su
% polinomi di grado 5 o inferiore, quindi ha ordine 6.
% Questa è anche la formula di quadratura predefinita di FreeFEM.

    if nargin < 6
        i = 1:cells.nc;
    end
    
    % Le coordinate dei punti di quadratura sono relative al triangolo
    % standard di vertici (0,0), (1,0), (0,1).
    qx = [1/3,      (6-sqrt(15))/21,     (6-sqrt(15))/21,     (9+2*sqrt(15))/21,...
                    (6+sqrt(15))/21,     (6+sqrt(15))/21,     (9-2*sqrt(15))/21];
    qy = [1/3,      (6-sqrt(15))/21,     (9+2*sqrt(15))/21,   (6-sqrt(15))/21,...
                    (6+sqrt(15))/21,     (9-2*sqrt(15))/21,   (6+sqrt(15))/21];
    qw = [270/1200, (155-sqrt(15))/1200, (155-sqrt(15))/1200, (155-sqrt(15))/1200,...
                    (155+sqrt(15))/1200, (155+sqrt(15))/1200, (155+sqrt(15))/1200];
    
    q_polygon = zeros(length(i),ndim);
    for j = 1:cells.mne
        e = abs(cells.e(i,j));
        mask = (e ~= 0);
        e = e(mask);
        v1 = edges.v1(e);
        v2 = edges.v2(e);
        ax = vertices.x(v1);
        ay = vertices.y(v1);
        bx = vertices.x(v2);
        by = vertices.y(v2);
        cx = cells.cx(i(mask));
        cy = cells.cy(i(mask));
        bax = bx-ax;
        bay = by-ay;
        cax = cx-ax;
        cay = cy-ay;
        area_triangle = 0.5*abs(bax.*cay - bay.*cax);
        q_triangle = zeros(sum(mask),ndim);
        for k = 1:length(qw)
            x = bax*qx(k) + cax*qy(k) + ax;
            y = bay*qx(k) + cay*qy(k) + ay;
            q_triangle = q_triangle + qw(k)*g(x,y);
        end
        q_triangle = area_triangle .* q_triangle;
        q_polygon(mask,:) = q_polygon(mask,:) + q_triangle;
    end
end



