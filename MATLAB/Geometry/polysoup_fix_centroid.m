function [polysoup] = polysoup_fix_centroid(polysoup)
%POLYSOUP_FIX_CENTROID Ricalcola il baricentro di ogni poligono.
    polysoup.cx = polysoup.vx(polysoup.p(:,1));
    polysoup.cy = polysoup.vy(polysoup.p(:,1));
    g = @(x,y) 1;
    gx = @(x,y) x;
    gy = @(x,y) y;
    area = polysoup_integral(g,1,polysoup);
    cx = polysoup_integral(gx,1,polysoup)./area;
    cy = polysoup_integral(gy,1,polysoup)./area;
    polysoup.cx = cx;
    polysoup.cy = cy;
end
