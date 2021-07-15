function [p_out] = polysoup_clip(p_in,clipping_region,mnv_limit)
%POLYSOUP_CLIP Interseca la polysoup p_in con la polyshape clipping_region.
% La clipping region dev'essere interamente contenuta nella polysoup.
    
    % ===== PSEUDOCODICE =====
    % per ogni boundary della clipping region:
    %   scegli un segmento a caso
    %   trova il poligono p_i di p con cui si interseca (forza bruta)
    %   marca tutti i poligoni lungo questo bordo a partire da p_i
    % per ogni poligono di p:
    %   se non è marcato ed è esterno alla clipping region:
    %       flood fill che marca i poligoni come esterni
    %   se non è marcato ed è interno alla clipping region:
    %       flood fill che marca i poligoni come interni
    % per ogni poligono di p:
    %   se è interno:
    %       aggiungilo così com'è a p_clipped
    %   se è di bordo:
    %       ritaglialo rispetto alla clipping_region
    %   se è esterno:
    %       ignoralo
    
    if nargin < 3
        mnv_limit = 16;
    end
    
    % Marca tutti poligoni che intersecano il bordo della clipping_region
    % 0 -> non marcato
    % 1 -> interseca il bordo della clipping_region
    % 2 -> è interno alla clipping region
    % 3 -> è esterno alla clipping region
    [vertices,edges,cells] = polymesh_from_polysoup(p_in);
    marker = zeros(cells.nc,1,'uint8');
    for b = 1:clipping_region.numboundaries
        [crx,cry] = clipping_region.boundary(b);
        marker = max(marker,...
            polysoup_clip_mark_boundary(vertices,edges,cells,crx,cry));
    end
    
    % Marca i poligoni esterni e interni
    for i = 1:cells.nc
        if marker(i) ~= 0
            continue;
        end
        if clipping_region.isinterior(cells.cx(i),cells.cy(i))
            marker(i) = 2;
            marker = polysoup_clip_flood_fill(vertices,edges,cells,marker,i);
        else
            marker(i) = 3;
            marker = polysoup_clip_flood_fill(vertices,edges,cells,marker,i);
        end
    end
    
    % Raggruppa in una polyshape tutti i poligoni con marker = 1
    idc = find(marker == 1)';
    nvc = 0;
    npc = length(idc);
    for i = idc
        nvc = nvc + cells.ne(i);
    end
    psx = zeros(nvc+npc,1);
    psy = zeros(nvc+npc,1);
    counter = 1;
    for i = idc
        for j = 1:cells.ne(i)
            e = cells.e(i,j);
            v = 0;
            if e > 0
                v = edges.v1(e);
            elseif e < 0
                v = edges.v2(-e);
            end
            psx(counter) = vertices.x(v);
            psy(counter) = vertices.y(v);
            counter = counter + 1;
        end
        psx(counter) = NaN;
        psy(counter) = NaN;
        counter = counter + 1;
    end
    assert(nvc+npc == counter - 1);
    ps = polyshape(psx,psy,'Simplify',false);
    
    % Interseca ps con la clipping region
    psa = intersect(regions(ps),clipping_region);
    psa = polyshape_flatten_array(psa);
    
    % Calcolo del numero di vertici e poligoni in psa
    psa_np = length(psa);
    psa_nvp = zeros(psa_np,1);
    for i = 1:psa_np
        psa_nvp(i) = size(psa(i).Vertices,1);
    end
    psa_nv = sum(psa_nvp);
    psa_mnv = max(psa_nvp);
    
    % Calcolo del numero di vertici e poligoni in p_out
    nv = psa_nv;
    np = psa_np;
    mnv = psa_mnv;
    for i = find(marker == 2)'
        nv = nv + cells.ne(i);
        np = np + 1;
        mnv = max(mnv, cells.ne(i));
    end
    if mnv > mnv_limit
        error(['Il processo di clipping ha superato mnv_limit.', ...
            ' Riprova con mnv_limit = %u.'], mnv);
    end
    
    % Alloca la memoria per p_out
    p_out = struct();
    p_out.nv = nv;
    p_out.vx = zeros(nv,1);
    p_out.vy = zeros(nv,1);
    p_out.np = np;
    p_out.mnv = double(mnv);
    p_out.p = zeros(np,mnv,'uint32');
    p_out.cx = zeros(np,1);
    p_out.cy = zeros(np,1);
    
    % Popola p_out con i poligoni interni che non sono stati ritagliati
    counter_v = 1;
    counter_p = 1;
    for i = 1:cells.nc
        if marker(i) == 2
            for j = 1:cells.ne(i)
                e = cells.e(i,j);
                v = 0;
                if e > 0
                    v = edges.v1(e);
                elseif e < 0
                    v = edges.v2(-e);
                end
                p_out.vx(counter_v) = vertices.x(v);
                p_out.vy(counter_v) = vertices.y(v);
                p_out.p(counter_p,j) = counter_v;
                counter_v = counter_v + 1;
            end
            counter_p = counter_p + 1;
        elseif marker(i) == 0
            error('Alcune celle non sono state marcate');
        end
    end
    
    % Popola p_out con i poligoni che sono stati ritagliati (quelli di psa)
    for i = 1:psa_np
        psxy = psa(i).Vertices;
        l = size(psxy,1);
        range = counter_v : counter_v+l-1;
        p_out.vx(range) = psxy(:,1);
        p_out.vy(range) = psxy(:,2);
        p_out.p(counter_p,1:l) = range;
        counter_v = counter_v + l;
        counter_p = counter_p + 1;
    end
    assert(nv == counter_v - 1);
    assert(np == counter_p - 1);
    
    % Ritocchi finali    
    p_out = polysoup_fix_CCW(p_out);
    p_out = polysoup_fix_centroid(p_out);
end
