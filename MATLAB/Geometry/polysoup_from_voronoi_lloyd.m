function [polysoup] = polysoup_from_voronoi_lloyd(x,y,clipping_region,k)
%POLYSOUP_FROM_VORONOI_LLOYD Come la funzione polysoup_from_voronoi(),
% ma il diagramma di Voronoi viene ricostruito k volte a partire dal baricentro
% di ogni poligono (algoritmo di Lloyd). Questo consente di aggiungere al
% vettore dei nodi anche i vertici della clipping_region, così il bordo viene
% preservato in modo migliore.
    
    assert(size(x,2)==1);
    assert(size(y,2)==1);

    % Scarta i nodi in ingresso al di fuori della clipping region
    mask = isinterior(clipping_region,x,y);
    x = x(mask);
    y = y(mask);
    
    % Aggiunge i vertici della clipping_region ai nodi in ingresso
    vx = clipping_region.Vertices(:,1);
    vy = clipping_region.Vertices(:,2);
    mask = ~isnan(vx);
    x = [x; vx(mask)];
    y = [y; vy(mask)];
    C = uniquetol([x,y],'ByRows',true,'DataScale',1);
    x = C(:,1);
    y = C(:,2);
    
    % Calcola la polygonal soup a partire dal diagramma di Voronoi, poi
    % esegui le prime k-1 iterazioni dell'algoritmo di Lloyd con una limitazione
    % larga su mnv e l'ultima iterazione con una limitazione stretta.
    if k == 0
        polysoup = polysoup_from_voronoi(x,y,clipping_region);
    else
        polysoup = polysoup_from_voronoi(x,y,clipping_region,128);
        for i = 1:k-1
            polysoup = polysoup_from_voronoi(polysoup.cx,polysoup.cy,clipping_region,128);
        end
        polysoup = polysoup_from_voronoi(polysoup.cx,polysoup.cy,clipping_region);
    end
end






