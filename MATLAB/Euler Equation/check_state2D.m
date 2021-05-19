function [] = check_state2D(u)
%CHECK_STATE2D Controlla che lo stato u del fluido sia ammissibile.
% Il controllo di densità e pressione assicura che anche l'energia e la
% temperatura siano positive. Sulla velocità non ci sono limitazioni.
    if any(u(:,1) <= 0)
        error('Si è formato il vuoto.');
    end
    p = pressure2D(u);
    if any(p <= 0)
        error('Si è formata una zona a pressione nulla o negativa.');
    end
end
