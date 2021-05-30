function [] = vtk_from_polymesh_ascii(filename,vertices,edges,cells)
%VTK_FROM_POLYMESH_ASCII Salva la soluzione in un file vtk in formato testuale
    fileID = fopen(filename,'w');
    fprintf(fileID,'# vtk DataFile Version 4.2\n');
    fprintf(fileID,'vtk output\n');
    fprintf(fileID,'ASCII\n');
    fprintf(fileID,'DATASET POLYDATA\n');
    
    % Vertici
    fprintf(fileID,'POINTS %u double\n', uint32(vertices.nv));
    fprintf(fileID,'%.12e %.12e 0.0e+00\n', [vertices.x, vertices.y]');
    
    % Celle
    v = zeros(cells.nc,cells.mne,'uint32');
    for j = 1:cells.mne
        e = cells.e(:,j);
        maskp = (e>0);
        maskm = (e<0);
        v(maskp,j) = edges.v1(e(maskp));
        v(maskm,j) = edges.v2(-e(maskm));
    end
    fprintf(fileID,'POLYGONS %u %u\n', uint32(cells.nc), uint32(sum(1+cells.ne)));
    for i=1:cells.nc
        % Importante: gli indici partono da 0 nei file vtk
        fprintf(fileID,'%u ',uint32(cells.ne(i)),v(i,1:cells.ne(i))-1);
        fprintf(fileID,'\n');
    end
    
    fprintf(fileID,'CELL_DATA %u\n', uint32(cells.nc));
    
    % Densità di massa
    fprintf(fileID,'SCALARS MassDensity double 1\n');
    fprintf(fileID,'LOOKUP_TABLE default\n');
    fprintf(fileID,'%.12e\n',cells.u(:,1));
    
    % Densità della quantità di moto
    fprintf(fileID,'VECTORS MomentumDensity double\n');
    fprintf(fileID,'%.12e %.12e 0.0e+00\n',[cells.u(:,2),cells.u(:,3)]');
    
    % Densità dell'energia totale
    fprintf(fileID,'SCALARS TotalEnergyDensity double 1\n');
    fprintf(fileID,'LOOKUP_TABLE default\n');
    fprintf(fileID,'%.12e\n',cells.u(:,4));
    
    % Pressione
    fprintf(fileID,'SCALARS Pressure double 1\n');
    fprintf(fileID,'LOOKUP_TABLE default\n');
    fprintf(fileID,'%.12e\n',pressure2D(cells.u));
    
    % Velocità
    fprintf(fileID,'VECTORS Velocity double\n');
    vx = cells.u(:,2)./cells.u(:,1);
    vy = cells.u(:,3)./cells.u(:,1);
    fprintf(fileID,'%.12e %.12e 0.0e+00\n',[vx,vy]');
    
    % Temperatura
    fprintf(fileID,'SCALARS Temperature double 1\n');
    fprintf(fileID,'LOOKUP_TABLE default\n');
    fprintf(fileID,'%.12e\n',temperature2D(cells.u));
    
    % Numero di Mach
    fprintf(fileID,'SCALARS MachNumber double 1\n');
    fprintf(fileID,'LOOKUP_TABLE default\n');
    fprintf(fileID,'%.12e\n',mach_number2D(cells.u));
    
    fclose(fileID);
end







