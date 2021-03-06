if isunix
    delete *.mexw64
else
    delete *.mexa64
end

mex -R2018a CXXFLAGS='$CXXFLAGS -Wall' ...
    CXXLIBS='$CXXLIBS' ...
    CXXOPTIMFLAGS='$CXXOPTIMFLAGS -O3 -march=native' ...
    polymesh_dfb.cpp

mex -R2018a CXXFLAGS='$CXXFLAGS -Wall' ...
    CXXLIBS='$CXXLIBS' ...
    CXXOPTIMFLAGS='$CXXOPTIMFLAGS -O3 -march=native' ...
    polymesh_from_polysoup_merge_edges.cpp

mex -R2018a CXXFLAGS='$CXXFLAGS -Wall' ...
    CXXLIBS='$CXXLIBS' ...
    CXXOPTIMFLAGS='$CXXOPTIMFLAGS -O3 -march=native' ...
    polysoup_merge_vertices_cleanup.cpp

mex -R2018a CXXFLAGS='$CXXFLAGS -Wall' ...
    CXXLIBS='$CXXLIBS' ...
    CXXOPTIMFLAGS='$CXXOPTIMFLAGS -O3 -march=native' ...
    polysoup_clip_mark_boundary.cpp

mex -R2018a CXXFLAGS='$CXXFLAGS -Wall' ...
    CXXLIBS='$CXXLIBS' ...
    CXXOPTIMFLAGS='$CXXOPTIMFLAGS -O3 -march=native' ...
    polysoup_clip_flood_fill.cpp
