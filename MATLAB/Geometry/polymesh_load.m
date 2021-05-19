function [vertices,edges,cells] = polymesh_load(filename)
%POLYMESH_LOAD Carica una polygonal mesh da un file .mat
    load(filename,'vertices','edges','cells');
end
