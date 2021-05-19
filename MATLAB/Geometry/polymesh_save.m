function [] = polymesh_save(filename,vertices,edges,cells)
%POLYMESH_SAVE Salva una polygonal mesh in un file .mat
    save(filename,'vertices','edges','cells');
end
