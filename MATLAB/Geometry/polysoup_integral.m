function [q_polygon] = polysoup_integral(g,ndim,polysoup,i)
%POLYSOUP_INTEGRAL Calcola l'integrale di una funzione vettoriale g(x,y)
% di ndim componenti sui poligoni con indici i (dove "i" può essere un vettore).
% Ricordiamo che i poligoni devono essere necessariamente convessi
% e devono contenere sul bordo o al loro interno il punto (cells.cx,cells.cy).
% Il risultato q_polygon è il valore dell'integrale su ogni poligono;
% questi integrali locali possono poi essere sommati per ottenere
% un integrale globale.
% L'integrale su ogni poligono viene calcolato tramite
% una triangolazione "a ventaglio" intorno al centro,
% seguita da formule di quadratura su ogni triangolo del ventaglio.
% La formula per triangoli, dovuta a Radon, usa 7 punti ed è esatta su
% polinomi di grado 5 o inferiore, quindi ha ordine 6.
% Questa è anche la formula di quadratura predefinita di FreeFEM.

    if nargin < 4
        i = 1:polysoup.np;
    end
    
    % Le coordinate dei punti di quadratura sono relative al triangolo
    % standard di vertici (0,0), (1,0), (0,1).
    qx = [1/3,      (6-sqrt(15))/21,     (6-sqrt(15))/21,     (9+2*sqrt(15))/21,...
                    (6+sqrt(15))/21,     (6+sqrt(15))/21,     (9-2*sqrt(15))/21];
    qy = [1/3,      (6-sqrt(15))/21,     (9+2*sqrt(15))/21,   (6-sqrt(15))/21,...
                    (6+sqrt(15))/21,     (9-2*sqrt(15))/21,   (6+sqrt(15))/21];
    qw = [270/1200, (155-sqrt(15))/1200, (155-sqrt(15))/1200, (155-sqrt(15))/1200,...
                    (155+sqrt(15))/1200, (155+sqrt(15))/1200, (155+sqrt(15))/1200];
    
    nvp = sum(polysoup.p~=0,2)';
    q_polygon = zeros(length(i),ndim);
    for j = 1:polysoup.mnv
        v1 = polysoup.p(i,j);
        mask = (v1 ~= 0);
        v1 = v1(mask);
        jnext = mod(j,nvp(i(mask)))+1;
        v2 = polysoup.p(sub2ind(size(polysoup.p),i(mask),jnext));
        ax = polysoup.vx(v1);
        ay = polysoup.vy(v1);
        bx = polysoup.vx(v2);
        by = polysoup.vy(v2);
        cx = polysoup.cx(i(mask));
        cy = polysoup.cy(i(mask));
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

