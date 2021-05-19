function [p] = pressure2D(u)
    %PRESSURE2D Calcola la pressione di un gas nello stato u.
    p = (adiabatic_index-1)*( u(:,4) - 0.5*(u(:,2).*u(:,2)+u(:,3).*u(:,3))./u(:,1) );
end
