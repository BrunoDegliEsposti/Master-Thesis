function [ps] = polyshape_from_polysoup(polysoup)
%POLYSHAPE_FROM_POLYSOUP Costruisci una polyshape a partire da una polygonal soup.
    X = cell(1,polysoup.np);
    Y = cell(1,polysoup.np);
    for i=1:polysoup.np
        p = nonzeros(polysoup.p(i,:));
        X{i} = polysoup.vx(p);
        Y{i} = polysoup.vy(p);
    end
    ps = polyshape(X,Y,'Simplify',false);
    ps = simplify(ps,'KeepCollinearPoints',false);
end
