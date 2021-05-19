function [polysoup] = polysoup_from_voronoi_convdec(x,y,clipping_region)
%POLYSOUP_FROM_VORONOI_CONVDEC Costruisci una polysoup associata a un diagramma
% di Voronoi. Le celle del diagramma sono intersecate con la polyshape clipping_region.
% Le componenti connesse che si ottengono da ogni cella dopo l'intersezione
% vengono decomposte in regioni convesse tramite triangolazione.

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
    assert(counter==nv+np);
    [psx,psy] = poly2ccw(psx,psy);
    
    % Crea un array di polyshapes, interseca ogni elemento dell'array
    % con la clipping region e scarta le polyshapes vuote.
    ps = polyshape(psx,psy,'SolidBoundaryOrientation','ccw','Simplify',false);
    assert(all(size(clipping_region)==[1,1]));
    psa = intersect(regions(ps),clipping_region);
    mask = ones(size(psa),'logical');
    for i = 1:length(psa)
        ps = psa(i);
        if ps.NumRegions == 0
            mask(i) = false;
        end
    end
    psa = psa(mask);
    
    % Stima del numero di vertici e poligoni che andremo a creare.
    % La stima è in eccesso, supponendo il caso pessimo in cui triangoliamo tutto.
    estimate_nv = 0;
    estimate_np = 0;
    for i = 1:length(psa)
        ps = psa(i);
        estimate_nv = estimate_nv + 3*(size(ps.Vertices,1)-2);
        estimate_np = estimate_np + size(ps.Vertices,1)-2;
    end
    
    % Crea i dati per la polysoup
    vx = zeros(estimate_nv,1);
    vy = zeros(estimate_nv,1);
    nvp = zeros(estimate_np,1);
    vcounter = 0;
    pcounter = 0;
    for i = 1:length(psa)
        psa_local = regions(psa(i));
        for j = 1:length(psa_local)
            ps = psa_local(j);
            psvx = ps.Vertices(:,1);
            psvy = ps.Vertices(:,2);
            [psvx,psvy] = poly2ccw(psvx,psvy);
            psnv = length(psvx);
            tricky = any(isnan(psvx));
            tricky = tricky || (psnv>16);
            tricky = tricky || ~ispolyconvex(psvx,psvy);
            if tricky
                TR = triangulation(ps);
                cl = TR.ConnectivityList;
                for k = 1:size(cl,1)
                    vx_local = TR.Points(cl(k,:),1);
                    vy_local = TR.Points(cl(k,:),2);
                    if ispolycw(vx_local,vy_local)
                        vx(vcounter+[3,2,1]) = vx_local;
                        vy(vcounter+[3,2,1]) = vy_local;
                    else
                        vx(vcounter+[1,2,3]) = vx_local;
                        vy(vcounter+[1,2,3]) = vy_local;
                    end
                    nvp(pcounter+1) = 3;
                    vcounter = vcounter + 3;
                    pcounter = pcounter + 1;
                end
            else
                vx(vcounter+1:vcounter+psnv) = psvx;
                vy(vcounter+1:vcounter+psnv) = psvy;
                nvp(pcounter+1) = psnv;
                vcounter = vcounter + psnv;
                pcounter = pcounter + 1;
            end
        end
    end
    vx = vx(1:vcounter);
    vy = vy(1:vcounter);
    nvp = nvp(1:pcounter);
    
    % Costruisci la polysoup    
    polysoup = struct();
    polysoup.nv = vcounter;
    polysoup.vx = vx;
    polysoup.vy = vy;
    polysoup.np = pcounter;
    polysoup.mnv = max(nvp);
    polysoup.p = zeros(pcounter,max(nvp),'uint32');
    k = 1;
    for i = 1:pcounter
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






