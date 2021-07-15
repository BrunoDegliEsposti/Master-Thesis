function [p] = polysoup_from_voronoi_variable(x,y,clipping_region,k,mnv_limit)
%POLYSOUP_FROM_VORONOI_VARIABLE Costruisci una polysoup associata a un
% diagramma di Voronoi. Le celle del diagramma sono intersecate con la
% polyshape clipping_region. Il diagramma di Voronoi è regolarizzato k volte
% mediante l'agoritmo di Lloyd. Se il valore in uscita polysoup.mnv è
% maggiore di mnv_limit, il processo di costruzione della mesh viene
% interrotto. A ogni iterazione dell'algoritmo di Lloyd, nuovi vertici
% vengono introdotti lungo il bordo della clipping region.

    assert(size(x,2)==1);
    assert(size(y,2)==1);

    if nargin < 4
        k = 0;
    end
    assert(k >= 0 && k == round(k));

    if nargin < 5
        mnv_limit = 16;
    end
    
    [x,y] = add_cr_vertices(x,y,clipping_region);
    p = polysoup_from_voronoi_nodes(x,y,128);
    
    for i = 1:k
        p = polysoup_clip(p,clipping_region,128);
        if i < round(k*0.7)
            [x,y] = add_cr_vertices(p.cx,p.cy,clipping_region);
        else
            x = p.cx;
            y = p.cy;
        end
        p = polysoup_from_voronoi_nodes(x,y,128);
    end
    
    p = polysoup_clip(p,clipping_region,mnv_limit);
end



function [x,y] = add_cr_vertices(x,y,clipping_region)
    crx = clipping_region.Vertices(:,1);
    cry = clipping_region.Vertices(:,2);
    mask = ~isnan(crx);
    x = [x; crx(mask)];
    y = [y; cry(mask)];
    C = uniquetol([x,y],'ByRows',true,'DataScale',1);
    x = C(:,1);
    y = C(:,2);
end








