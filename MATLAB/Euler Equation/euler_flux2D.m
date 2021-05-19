function [f] = euler_flux2D(u)
    %EULER_FLUX2D Calcola il flusso dell'equazione di Eulero 2D (matrice 4x2).
    [n,~] = size(u);
    f = zeros(n,4,2);
    u1 = u(:,1);
    u2 = u(:,2);
    u3 = u(:,3);
    u4 = u(:,4);
    p = pressure2D(u);
    f(:,1,1) = u2;              f(:,1,2) = u3;
    f(:,2,1) = u2.*u2./u1 + p;  f(:,2,2) = u2.*u3./u1;
    f(:,3,1) = u3.*u2./u1;      f(:,3,2) = u3.*u3./u1 + p;
    f(:,4,1) = (u4+p).*u2./u1;  f(:,4,2) = (u4+p).*u3./u1;
end
