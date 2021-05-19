%% Regular disk
addpath('../Geometry');
cx = 0;
cy = 0;
radius = 1;

Nr = 20;
p = polysoup_from_grid(round(2*pi*Nr),Nr,0,0,2*pi,radius);
p = polysoup_transform(p,@pol2cart);
p = polysoup_transform(p,@(x,y) deal(x+cx,y+cy));
[v,e,c] = polymesh_from_polysoup(p);
polymesh_save('regular_disk_20.mat',v,e,c);

Nr = 40;
p = polysoup_from_grid(round(2*pi*Nr),Nr,0,0,2*pi,radius);
p = polysoup_transform(p,@pol2cart);
p = polysoup_transform(p,@(x,y) deal(x+cx,y+cy));
[v,e,c] = polymesh_from_polysoup(p);
polymesh_save('regular_disk_40.mat',v,e,c);

Nr = 80;
p = polysoup_from_grid(round(2*pi*Nr),Nr,0,0,2*pi,radius);
p = polysoup_transform(p,@pol2cart);
p = polysoup_transform(p,@(x,y) deal(x+cx,y+cy));
[v,e,c] = polymesh_from_polysoup(p);
polymesh_save('regular_disk_80.mat',v,e,c);



