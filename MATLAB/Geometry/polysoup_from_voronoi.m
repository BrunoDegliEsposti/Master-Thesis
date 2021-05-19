function [polysoup] = polysoup_from_voronoi(x,y,clipping_region,mnv_limit)
%POLYSOUP_FROM_VORONOI Costruisci una polysoup associata a un diagramma di Voronoi.
% Le celle del diagramma sono intersecate con la polyshape clipping_region.
% Le componenti connesse che si ottengono da ogni cella dopo l'intersezione
% vengono sostituite dal loro inviluppo convesso. Se il valore in uscita
% polysoup.mnv dovesse essere maggiore di mnv_limit, il processo di costruzione
% della mesh viene interrotto.

    if nargin < 4
        mnv_limit = 16;
    end

    % Costruisci il diagramma di Voronoi in modo che tutte le celle
    % associate ai nodi in ingresso siano limitate.
    r = 10*max(hypot(x,y));
    x(end+(1:8)) = [ r, r, r, 0,-r,-r,-r, 0];
    y(end+(1:8)) = [-r, 0, r, r, r, 0,-r,-r];
    DT = delaunayTriangulation(x,y);
    [v,c] = voronoiDiagram(DT);
    [np,~] = size(c);
    
    % Rimuovi celle vuote o illimitate
    mask = ones(np,1,'logical');
    for i = 1:np
        p = c{i};
        if (length(p) < 3) || any(p == 1)
            mask(i) = false;
        end
    end
    c = c(mask);
    np = sum(mask);
    
    % Crea un vettore di coordinate nel formato richiesto da polyshape()
    nv = 0;
    for i = 1:np
        nv = nv + length(c{i});
    end
    psx = zeros(nv+np,1);
    psy = zeros(nv+np,1);
    counter = 0;
    for i = 1:np
        p = c{i};
        l = length(p);
        psx(counter+1:counter+l) = v(p,1);
        psx(counter+l+1) = NaN;
        psy(counter+1:counter+l) = v(p,2);
        psy(counter+l+1) = NaN;
        counter = counter + l + 1;
    end
    
    % Crea un array di polyshapes, interseca ogni suo elemento con
    % la clipping region, scarta ogni poligono vuoto e infine
    % costruisci l'inviluppo convesso dei poligoni rimasti.
    ps = polyshape(psx,psy,'Simplify',false);
    assert(all(size(clipping_region)==[1,1]));
    psa = intersect(regions(ps),clipping_region);
    psa = polyshape_flatten_array(psa);
    psa = convhull(psa);
    
    % Calcolo del numero di poligoni e vertici che andremo a creare.
    % Grazie all'inviluppo convesso, la stima è esatta.
    np = length(psa);
    nvp = zeros(np,1);
    for i = 1:np
        ps = psa(i);
        assert(ps.NumRegions == 1);
        psnv = size(ps.Vertices,1);
        nvp(i) = psnv;
    end
    nv = sum(nvp);
    mnv = max(nvp);
    if mnv > mnv_limit
        error(['Almeno una cella contiene un numero eccessivo di vertici. ',...
            'Prova ad alzare mnv_limit a %u.'], mnv);
    end
    
    % Calcolo delle coordinate dei vertici della polysoup
    vx = zeros(nv,1);
    vy = zeros(nv,1);
    vcounter = 0;
    for i = 1:length(psa)
        ps = psa(i);
        psvx = ps.Vertices(:,1);
        psvy = ps.Vertices(:,2);
        [psvx,psvy] = poly2ccw(psvx,psvy);
        psnv = length(psvx);
        vx(vcounter+1:vcounter+psnv) = psvx;
        vy(vcounter+1:vcounter+psnv) = psvy;
        vcounter = vcounter + psnv;
    end
    assert(vcounter==nv);
    
    % Costruisci la polysoup    
    polysoup = struct();
    polysoup.nv = nv;
    polysoup.vx = vx;
    polysoup.vy = vy;
    polysoup.np = np;
    polysoup.mnv = mnv;
    polysoup.p = zeros(np,mnv,'uint32');
    k = 1;
    for i = 1:np
        for j = 1:nvp(i)
            polysoup.p(i,j) = k;
            k = k + 1;
        end
    end
    
    % Calcola il centro di messa di ogni cella
    polysoup.cx = zeros(np,1);
    polysoup.cy = zeros(np,1);
    polysoup = polysoup_fix_centers(polysoup);
end





