function [polysoup] = polysoup_from_voronoi_old(x,y,clipping_region,lloyd)
%POLYSOUP_FROM_VORONOI_OLD Costruisci una polysoup associata a un diagramma di Voronoi.
% Le celle del diagramma sono intersecate con la polyshape clipping_region.
% Lloyd è un argomento facoltativo di tipo intero e indica il numero di
% iterazioni dell'algoritmo di Lloyd da eseguire sul diagramma di Voronoi.
    
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
        l = length(p);
        if (l < 3) || any(p == 1)
            mask(i) = false;
        end
    end
    c = c(mask);
    np = sum(mask);
    
    % Interseca le celle con la clipping_region
    vx = zeros(16*np,1);
    vy = zeros(16*np,1);
    p = cell(np,1);
    cx = zeros(np,1);
    cy = zeros(np,1);
    vcounter = 1;
    pcounter = 1;
    mnv = 3;
    for i = 1:np
        % Intersezione e inviluppo convesso
        px = v(c{i},1);
        py = v(c{i},2);
        ps = polyshape(px,py,'Simplify',false);
        ps = intersect(ps,clipping_region);
        ps = convhull(ps);
        ps = simplify(ps,'KeepCollinearPoints',false);
        if ps.NumRegions == 0
            continue;
        end
        assert(ps.NumRegions == 1);
        % Costruzione del vettore dei vertici
        nv = size(ps.Vertices,1);
        assert(nv <= 16);
        mnv = max(mnv,nv);
        p{pcounter} = vcounter+(0:nv-1);
        for j = 1:nv
            vx(vcounter) = ps.Vertices(j,1);
            vy(vcounter) = ps.Vertices(j,2);
            vcounter = vcounter + 1;
        end
        % Ricostruzione del nodo che ha generato questo poligono di Voronoi
        [cenx,ceny] = centroid(ps);
        k = DT.nearestNeighbor(cenx,ceny);
        if isinterior(ps,x(k),y(k))
            cx(pcounter) = x(k);
            cy(pcounter) = y(k);
        else
            cx(pcounter) = cenx;
            cy(pcounter) = ceny;
        end
        pcounter = pcounter + 1;
    end
    nv = vcounter-1;
    vx = vx(1:nv);
    vy = vy(1:nv);
    np = pcounter-1;
    p = p(1:np);
    cx = cx(1:np);
    cy = cy(1:np);
    
    % Allocazione di memoria per la polysoup
    polysoup = struct();
    polysoup.nv = nv;
    polysoup.vx = vx;
    polysoup.vy = vy;
    polysoup.np = np;
    polysoup.mnv = mnv;
    polysoup.p = zeros(np,mnv,'uint32');
    polysoup.cx = cx;
    polysoup.cy = cy;
    
    % Converti il cell array p nella matrice polysoup.p
    for i = 1:np
        l = length(p{i});
        polysoup.p(i,1:l) = p{i};
    end
    
    % Rilassamento di Lloyd
    if (nargin > 3) && (lloyd > 0)
        polysoup = polysoup_fix_centers(polysoup);
        polysoup = polysoup_from_voronoi(polysoup.cx,polysoup.cy,clipping_region,lloyd-1);
    end
    
    % Aggiusta l'orientamento dei poligoni (se necessario)
    polysoup = polysoup_fix_CCW(polysoup);
end





