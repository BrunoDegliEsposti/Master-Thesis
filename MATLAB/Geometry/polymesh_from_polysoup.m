function [vertices,edges,cells] = polymesh_from_polysoup(polysoup,tol)
%POLYMESH_FROM_POLYSOUP Converte una polygonal soup in una mesh poligonale
% vera e propria, andando a dedurre la connettività mancante.

    if nargin < 2
        tol = 1e-8;
    end

    polysoup = polysoup_merge_vertices(polysoup,tol);
    
    vertices = struct();
    vertices.nv = polysoup.nv;
    vertices.x = polysoup.vx;
    vertices.y = polysoup.vy;
    cells = struct();
    cells.nc = polysoup.np;
    cells.cx = polysoup.cx;
    cells.cy = polysoup.cy;
    cells.mne = polysoup.mnv;
    cells.ne = zeros(polysoup.np,1,'uint8');
    cells.e = zeros(size(polysoup.p),'int32');
    cells.nac = zeros(polysoup.np,1,'uint8');
    cells.ac = zeros(size(polysoup.p),'uint32');
    edges = struct();
    edges.v1 = zeros(polysoup.np*polysoup.mnv,1,'uint32');
    edges.v2 = zeros(polysoup.np*polysoup.mnv,1,'uint32');
    edges.type = zeros(polysoup.np*polysoup.mnv,1,'uint32');
    edges.cp = zeros(polysoup.np*polysoup.mnv,1,'uint32');
    edges.cm = zeros(polysoup.np*polysoup.mnv,1,'uint32');
    
    % Ricostruisci gli spigoli
    hash_table = containers.Map('KeyType','uint64','ValueType','uint32');
    ecounter = 1;
    for i = 1:polysoup.np
        ne = sum(polysoup.p(i,:)~=0);
        cells.ne(i) = ne;
        for j = 1:ne
            jnext = mod(j,ne)+1;
            v1 = polysoup.p(i,j);
            v2 = polysoup.p(i,jnext);
            v12 = bitmerge(min(v1,v2),max(v1,v2));
            if isKey(hash_table,v12)
                k = hash_table(v12);
                edges.type(k) = edges.type(k) + 1;
                cells.e(i,j) = -int32(k);
                edges.cm(k) = i;
            else
                hash_table(v12) = ecounter;
                edges.v1(ecounter) = v1;
                edges.v2(ecounter) = v2;
                edges.type(ecounter) = edges.type(ecounter) + 1;
                cells.e(i,j) = ecounter;
                edges.cp(ecounter) = i;
                ecounter = ecounter + 1;
            end
        end
    end
    edges.ne = ecounter-1;
    edges.v1 = edges.v1(1:edges.ne);
    edges.v2 = edges.v2(1:edges.ne);
    edges.type = edges.type(1:edges.ne);
    edges.cp = edges.cp(1:edges.ne);
    edges.cm = edges.cm(1:edges.ne);
    if any(edges.type == 0) || any(edges.type > 2)
        error('Lo stesso spigolo è condiviso da 3 o più facce')
    end
    edges.type = 2-edges.type;
    edges.nbe = sum(edges.type);
    edges.nie = edges.ne - edges.nbe;
    edges.length = edge_length(vertices,edges);
    assert(all(edges.length>0));
    [edges.nx,edges.ny] = edge_normal(vertices,edges);
    
    % Ricostruisci l'adiacenza tra celle
    for i = 1:cells.nc
        for j = 1:cells.ne(i)
            k = cells.e(i,j);
            if k > 0
                cells.ac(i,j) = edges.cm(k);
            elseif k < 0
                cells.ac(i,j) = edges.cp(-k);
            end
        end
    end
    cells.nac = sum(cells.ac~=0,2);
    cells.area = cell_area(vertices,edges,cells);
    assert(all(cells.area>0));
    cells.perimeter = cell_perimeter(vertices,edges,cells);
    assert(all(cells.perimeter>0));
    
    % Riordina gli spigoli: prima quelli interni, poi quelli sul bordo
    P = [find(edges.type == 0); find(edges.type == 1)];
    edges.type = edges.type(P);
    edges.v1 = edges.v1(P);
    edges.v2 = edges.v2(P);
    edges.length = edges.length(P);
    edges.nx = edges.nx(P);
    edges.ny = edges.ny(P);
    edges.cp = edges.cp(P);
    edges.cm = edges.cm(P);
    invP(P) = 1:edges.ne;
    for i = 1:cells.nc
        for j = 1:cells.ne(i)
            k = cells.e(i,j);
            if k ~= 0
                cells.e(i,j) = sign(k)*invP(abs(k));
            end
        end
    end
end

function [u] = bitmerge(in1,in2)
    ul = bitshift(uint64(in1),32);
    ur = uint64(in2);
    u = bitor(ul,ur);
end



