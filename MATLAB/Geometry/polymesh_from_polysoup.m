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
    
    [edges, cells] = polymesh_from_polysoup_merge_edges(polysoup);
    cells.cx = polysoup.cx;
    cells.cy = polysoup.cy;
    
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
    cells.nac = zeros(cells.nc,1,'uint32');
    cells.ac = zeros(size(cells.e),'uint32');
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
    cells.nac = uint32(sum(cells.ac~=0,2));
    cells.area = cell_area(vertices,edges,cells);
    assert(all(cells.area>0));
    cells.perimeter = cell_perimeter(vertices,edges,cells);
    assert(all(cells.perimeter>0));
    cells.h = 4*cells.area./cells.perimeter;
    
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
    
    % Calcola la distanza di ogni cella dal bordo (numero di vicini)
    cells.dfb = polymesh_dfb(vertices,edges,cells);
end
