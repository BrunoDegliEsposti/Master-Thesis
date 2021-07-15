function [p] = polysoup_from_voronoi_uniform(x,y,clipping_region,k,mnv_limit)
%POLYSOUP_FROM_VORONOI_UNIFORM Costruisci una polysoup associata a un diagramma
% di Voronoi. Le celle del diagramma sono intersecate con la polyshape
% clipping_region. Il diagramma di Voronoi è regolarizzato k volte
% mediante l'agoritmo di Lloyd. Se il valore in uscita polysoup.mnv è
% maggiore di mnv_limit, il processo di costruzione della mesh viene
% interrotto. Le celle tendono a distribuirsi in modo uniforme su
% tutto il dominio.

    assert(size(x,2)==1);
    assert(size(y,2)==1);

    if nargin < 4
        k = 0;
    end
    assert(k >= 0 && k == round(k));

    if nargin < 5
        mnv_limit = 16;
    end
    
    p = polysoup_from_voronoi_nodes(x,y,128);
    
    for i = 1:k
        p = polysoup_clip(p,clipping_region,128);
        p = polysoup_from_voronoi_nodes(p.cx,p.cy,128);
    end
    
    p = polysoup_clip(p,clipping_region,mnv_limit);
end
