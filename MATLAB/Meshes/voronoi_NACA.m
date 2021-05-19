%% Profilo alare NACA 0012
addpath('../Geometry');
r = 5;
angle_of_attack = 4;
vxy = readmatrix('NACA0012-Selig.dat');
ps1 = polyshape(vxy(:,1),vxy(:,2));
ps1 = translate(ps1,-0.5,0);
ps1 = rotate(ps1,-angle_of_attack);

N = 2500;
hs = haltonset(2);
x = 2*r*hs(1:N,1)-r;
y = 2*r*hs(1:N,2)-r;
ps2 = nsidedpoly(round(sqrt(N)*pi),'Radius',r);
clipping_region = subtract(ps2,ps1);
% Nr = 20;
% Ntheta = 120;
% r = linspace(0,5,Nr);
% theta = linspace(0,2*pi,Ntheta);
% [R,THETA] = meshgrid(r,theta);
% [x,y] = pol2cart(THETA(:),R(:));
p = polysoup_from_voronoi_lloyd2(x,y,clipping_region,3);
[v,e,c] = polymesh_from_polysoup(p);
polymesh_save('voronoi_NACA_2500',v,e,c);

N = 10000;
hs = haltonset(2);
x = 2*r*hs(1:N,1)-r;
y = 2*r*hs(1:N,2)-r;
ps2 = nsidedpoly(round(sqrt(N)*pi),'Radius',r);
clipping_region = subtract(ps2,ps1);
% Nr = 20;
% Ntheta = 120;
% r = linspace(0,5,Nr);
% theta = linspace(0,2*pi,Ntheta);
% [R,THETA] = meshgrid(r,theta);
% [x,y] = pol2cart(THETA(:),R(:));
p = polysoup_from_voronoi_lloyd2(x,y,clipping_region,3);
[v,e,c] = polymesh_from_polysoup(p);
polymesh_save('voronoi_NACA_10000',v,e,c);


