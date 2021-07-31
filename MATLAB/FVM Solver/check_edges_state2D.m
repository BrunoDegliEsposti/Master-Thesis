function [] = check_edges_state2D(edges)
%CHECK_EDGES_STATE2D Controlla che gli stati um,up ricostruiti sui due lati
% di ogni spigolo siano ammissibili.
    for k = 1:edges.nq
        if any(edges.um(:,1,k) <= 0)
            error('Su uno spigolo è stata ricostruita una densità <= 0.');
        end
        if any(edges.up(:,1,k) <= 0)
            error('Su uno spigolo è stata ricostruita una densità <= 0.');
        end
        pm = pressure2D(edges.um(:,:,k));
        if any(pm <= 0)
            error('Su uno spigolo è stata ricostruita una pressione <= 0.');
        end
        pp = pressure2D(edges.up(:,:,k));
        if any(pp <= 0)
            error('Su uno spigolo è stata ricostruita una pressione <= 0.');
        end
    end
end

