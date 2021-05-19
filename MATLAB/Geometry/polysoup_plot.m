function [] = polysoup_plot(polysoup,i)
%POLYSOUP_PLOT Disegna una polygonal soup.
    if nargin == 1
        i = 1:polysoup.np;
    end
    polysoup.p = double(polysoup.p(i,:));
    mask = (polysoup.p == 0);
    polysoup.p(mask) = NaN;
    figure;
    hold on;
    patch('Faces',polysoup.p,'Vertices',[polysoup.vx,polysoup.vy],...
          'FaceColor','#DDD');
    scatter(polysoup.cx(i),polysoup.cy(i),'.k');
    hold off;
    axis equal;
end
