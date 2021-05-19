function [psa_flat] = polyshape_flatten_array(psa)
%POLYSHAPE_FLATTEN_ARRAY Converte un polyshape array in cui ogni elemento può
% contenere più regioni in un polyshape array in cui ogni elemento
% contiene una sola regione.

    % Alloca memoria
    nregions = [psa.NumRegions];
    len_flat = sum(nregions);
    psa_flat = polyshape.empty();
    psa_flat(len_flat,1) = polyshape();
    
    % Copia le regioni
    pscounter = 0;
    for i=1:length(psa)
        psa_local = regions(psa(i));
        for j = 1:length(psa_local)
            psa_flat(pscounter+1,1) = psa_local(j);
            pscounter = pscounter + 1;
        end
    end
    assert(pscounter==len_flat);
end
