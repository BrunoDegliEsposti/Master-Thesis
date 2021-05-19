function [polysoup] = polysoup_merge_vertices(polysoup,tol)
%POLYSOUP_MERGE_VERTICES Rimuove i vertici duplicati da una polygonal soup.
% Ogni cluster di vertici di raggio minore di tol viene accorpato
% in un unico vertice. Questo è un passo fondamentale per calcolare
% la connettività a partire da una polysoup.

    % Rimuovi i vertici duplicati
    [vxy_new,~,indices_new] = uniquetol(...
        [polysoup.vx,polysoup.vy],tol,'ByRows',true,'DataScale',1);
    polysoup.nv = size(vxy_new,1);
    polysoup.vx = vxy_new(:,1);
    polysoup.vy = vxy_new(:,2);
    
    % Aggiorna gli indici dei poligoni
    for j=1:polysoup.mnv
        v = polysoup.p(:,j);
        mask = (v~=0);
        v(mask) = indices_new(v(mask));
        polysoup.p(:,j) = v;
    end
    
    % Rimuovi i poligoni con meno di tre vertici distinti
    mnv = 3;
    mask = ones(polysoup.np,1,'logical');
    for i=1:polysoup.np
        v = polysoup.p(i,:);
        w = unique(v,'stable');
        lw = length(w);
        v = v*0;
        v(1:lw) = w;
        polysoup.p(i,:) = v;
        if lw < 3
            mask(i) = false;
        else
            mnv = max(mnv,lw);
        end
        % Volendo, qui si potrebbe controllare che il poligono sia
        % ancora convesso, che contenga (cx,cy), che i suoi vertici
        % non siano tutti allineati, ecc...
    end
    polysoup.np = sum(mask);
    polysoup.mnv = mnv;
    polysoup.p = polysoup.p(mask,1:mnv);
    polysoup.cx = polysoup.cx(mask);
    polysoup.cy = polysoup.cy(mask);
end



