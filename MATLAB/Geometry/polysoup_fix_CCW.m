function [polysoup] = polysoup_fix_CCW(polysoup)
%POLYSOUP_FIX_CCW Corregge l'orientazione di ogni poligono in una polygonal soup,
% affinché i vertici di ogni poligono siano elencati in ordine antiorario.
    ax = polysoup.vx(polysoup.p(:,1));
    ay = polysoup.vy(polysoup.p(:,1));
    bx = polysoup.vx(polysoup.p(:,2));
    by = polysoup.vy(polysoup.p(:,2));
    cx = polysoup.vx(polysoup.p(:,3));
    cy = polysoup.vy(polysoup.p(:,3));
    b = counterclockwise(ax,ay,bx,by,cx,cy);
    indices = find(not(b)');
    for i = indices
        l = sum(polysoup.p(i,:)~=0);
        polysoup.p(i,1:l) = fliplr(polysoup.p(i,1:l));
    end
end

