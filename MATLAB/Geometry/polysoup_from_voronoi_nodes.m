function [polysoup] = polysoup_from_voronoi_nodes(x,y,mnv_limit)
%POLYSOUP_FROM_VORONOI_NODES Costruisci una polysoup associata alle celle
% del diagramma di Voronoi generato dai nodi (x,y). La polysoup non
% è pronta per l'uso, ma va tagliata lungo il bordo.

    if nargin < 3
        mnv_limit = 16;
    end

    % Costruisci il diagramma di Voronoi in modo che tutte le celle
    % associate ai nodi in ingresso siano limitate
    r = 10*max(hypot(x,y));
    x(end+(1:8)) = [ r, r, r, 0,-r,-r,-r, 0];
    y(end+(1:8)) = [-r, 0, r, r, r, 0,-r,-r];
    [v,c] = voronoin([x,y]);
    nc = size(c,1);
    
    % Rileva celle vuote, illimitate o con troppi spigoli, affinché
    % si possano calcolare i valori esatti di nv, np e mnv.
    mask = ones(nc,1,'logical');
    nv = 0;
    np = 0;
    mnv = 0;
    for i = 1:nc
        p = c{i};
        l = length(p);
        if length(p) < 3
            % cella vuota
            mask(i) = false;
        elseif any(p == 1)
            % cella illimitata
            mask(i) = false;
        elseif l > mnv_limit
            % cella con troppi spigoli
            error("Una cella del diagramma di Voronoi ha troppi spigoli");
        else
            % cella valida
            nv = nv + l;
            np = np + 1;
            mnv = max(mnv, l);
        end
    end
    
    % Crea la polysoup, convertendo ogni cella valida c{i} in un poligono
    polysoup = struct();
    polysoup.nv = nv;
    polysoup.vx = zeros(nv,1);
    polysoup.vy = zeros(nv,1);
    polysoup.np = np;
    polysoup.mnv = mnv;
    polysoup.p = zeros(np,mnv,'uint32');
    counter_v = 1;
    counter_p = 1;
    for i = 1:nc
        if mask(i)
            p = c{i};
            l = length(p);
            polysoup.vx(counter_v:counter_v+l-1) = v(p,1);
            polysoup.vy(counter_v:counter_v+l-1) = v(p,2);
            polysoup.p(counter_p,1:l) = counter_v:counter_v+l-1;
            counter_v = counter_v + l;
            counter_p = counter_p + 1;
        end
    end
    assert(counter_v == nv+1);
    assert(counter_p == np+1);
    polysoup = polysoup_fix_CCW(polysoup);
    polysoup.cx = zeros(np,1);
    polysoup.cy = zeros(np,1);
    polysoup = polysoup_fix_centroid(polysoup);
end
