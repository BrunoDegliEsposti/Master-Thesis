%% Regular square
addpath('../Geometry');
blcx = -5;
blcy = -5;
Lx = 10;
Ly = 10;

Nx = 50;
Ny = 50;
p = polysoup_from_grid(Nx,Ny,blcx,blcy,Lx,Ly);
[v,e,c] = polymesh_from_polysoup(p);
polymesh_save('regular_square_50x50.mat',v,e,c);

Nx = 63;
Ny = 63;
p = polysoup_from_grid(Nx,Ny,blcx,blcy,Lx,Ly);
[v,e,c] = polymesh_from_polysoup(p);
polymesh_save('regular_square_63x63.mat',v,e,c);

Nx = 79;
Ny = 79;
p = polysoup_from_grid(Nx,Ny,blcx,blcy,Lx,Ly);
[v,e,c] = polymesh_from_polysoup(p);
polymesh_save('regular_square_79x79.mat',v,e,c);

Nx = 100;
Ny = 100;
p = polysoup_from_grid(Nx,Ny,blcx,blcy,Lx,Ly);
[v,e,c] = polymesh_from_polysoup(p);
polymesh_save('regular_square_100x100.mat',v,e,c);

Nx = 126;
Ny = 126;
p = polysoup_from_grid(Nx,Ny,blcx,blcy,Lx,Ly);
[v,e,c] = polymesh_from_polysoup(p);
polymesh_save('regular_square_126x126.mat',v,e,c);

Nx = 159;
Ny = 159;
p = polysoup_from_grid(Nx,Ny,blcx,blcy,Lx,Ly);
[v,e,c] = polymesh_from_polysoup(p);
polymesh_save('regular_square_159x159.mat',v,e,c);

Nx = 200;
Ny = 200;
p = polysoup_from_grid(Nx,Ny,blcx,blcy,Lx,Ly);
[v,e,c] = polymesh_from_polysoup(p);
polymesh_save('regular_square_200x200.mat',v,e,c);

Nx = 400;
Ny = 400;
p = polysoup_from_grid(Nx,Ny,blcx,blcy,Lx,Ly);
[v,e,c] = polymesh_from_polysoup(p);
polymesh_save('regular_square_400x400.mat',v,e,c);

Nx = 800;
Ny = 800;
p = polysoup_from_grid(Nx,Ny,blcx,blcy,Lx,Ly);
[v,e,c] = polymesh_from_polysoup(p);
polymesh_save('regular_square_800x800.mat',v,e,c);



