function [polysoup] = polysoup_fix_centers(polysoup)
%POLYSOUP_FIX_CENTERS Ricalcola i centri di ogni poligono come centro di massa.
    for i=1:polysoup.np
        p = nonzeros(polysoup.p(i,:));
        polysoup.cx(i) = mean(polysoup.vx(p));
        polysoup.cy(i) = mean(polysoup.vy(p));
    end
    g = @(x,y) 1;
    gx = @(x,y) x;
    gy = @(x,y) y;
    area = polysoup_integral(g,1,polysoup);
    cx = polysoup_integral(gx,1,polysoup)./area;
    cy = polysoup_integral(gy,1,polysoup)./area;
    polysoup.cx = cx;
    polysoup.cy = cy;
end
