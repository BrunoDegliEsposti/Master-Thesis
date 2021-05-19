function [] = polymesh_plot(vertices,edges,cells,i)
%POLYMESH_PLOT Disegna una polygonal mesh.
    if nargin < 4
        i = 1:cells.nc;
    elseif size(i,2) == 1
        i = i';
    end
    
    % Crea una matrice nel formato adatto a patch()
    v = zeros(length(i),cells.mne);
    v(:,:) = NaN;
    for j = 1:cells.mne
        e = cells.e(i,j);
        maskp = (e>0);
        maskm = (e<0);
        v(maskp,j) = edges.v1(e(maskp));
        v(maskm,j) = edges.v2(-e(maskm));
    end
    
    % Disegna le facce in griglio, gli spigoli in nero e i centri di ogni poligono
    figure;
    hold on;
    patch('Faces',v,'Vertices',[vertices.x,vertices.y],...
          'FaceColor','#DDD','DisplayName','Celle');
    scatter(cells.cx(i),cells.cy(i),'.k','DisplayName','Centri');
    
    % Disegna le normali a ogni spigolo
    j = cells.e(i,:);
    j = nonzeros(unique(abs(j(:))));
    [mx,my] = edge_midpoint(vertices,edges,j);
    [nx,ny] = edge_normal(vertices,edges,j);
    quiver(mx,my,nx,ny,0.2,'k','DisplayName','Normali');
    
    % Disegna gli spigoli sul bordo con colori diversi in base
    % al'identificativo delle condizioni al bordo
    j = edges.nie + (1:edges.nbe);
    types = unique(edges.type(j));
    for t = types'
        mask = (edges.type(j)==t);
        jt = j(mask);
        v1 = edges.v1(jt);
        v2 = edges.v2(jt);
        v1x = vertices.x(v1);
        v1y = vertices.y(v1);
        v2x = vertices.x(v2);
        v2y = vertices.y(v2);
        px = [v1x,v2x,NaN(size(v1x))]';
        py = [v1y,v2y,NaN(size(v1y))]';
        plot(px(:),py(:),'LineWidth',3,'DisplayName',['Bordo #',num2str(t)]);
    end
    
    legend;
    hold off;
    axis equal;
end




