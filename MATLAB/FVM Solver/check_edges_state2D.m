function [] = check_edges_state2D(edges)
%CHECK_EDGES_STATE2D Controlla che gli stati um,up ricostruiti sui due lati
% di ogni spigolo siano ammissibili.
    for k = 1:edges.nq
        if any(edges.um(:,1,k) <= 0)
            j = find(edges.um(:,1,k)<=0, 1);
            error('Su uno spigolo (j = %d, lato um) è stata ricostruita una densità <= 0.', j);
        end
        if any(edges.up(:,1,k) <= 0)
            j = find(edges.up(:,1,k)<=0, 1);
            error('Su uno spigolo (j = %d, lato up) è stata ricostruita una densità <= 0.', j);
        end
        pm = pressure2D(edges.um(:,:,k));
        if any(pm <= 0)
            j = find(pm<=0, 1);
            error('Su uno spigolo (j = %d, lato um) è stata ricostruita una pressione <= 0.', j);
        end
        pp = pressure2D(edges.up(:,:,k));
        if any(pp <= 0)
            j = find(pp<=0, 1);
            error('Su uno spigolo (j = %d, lato up) è stata ricostruita una pressione <= 0.', j);
        end
    end
end
