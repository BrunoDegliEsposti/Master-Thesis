% Funziona davvero, se raddoppio Nx e Ny l'errore scende di 2^6

addpath('../Geometry');
blcx = 0;
blcy = 0;
Lx = pi;
Ly = pi;
g = @(x,y) sin(x).*sin(y);

err = zeros(5,1);
for k=1:5
    Nx = 2^k;
    Ny = 2^k;
    p = polysoup_from_grid(Nx,Ny,blcx,blcy,Lx,Ly);
    q = polysoup_integral(g,1,p);
    err(k) = abs(sum(q)-4);
end

disp(err(1:4)./err(2:5));
