mex -R2018a CXXFLAGS='$CXXFLAGS -Wall' ...
    CXXLIBS='$CXXLIBS' ...
    CXXOPTIMFLAGS='$CXXOPTIMFLAGS -O3 -march=native' ...
    reconstruction_LLS1.cpp

mex -R2018a CXXFLAGS='$CXXFLAGS -Wall' ...
    CXXLIBS='$CXXLIBS' ...
    CXXOPTIMFLAGS='$CXXOPTIMFLAGS -O3 -march=native' ...
    reconstruction_LLS2_naive.cpp -lmwlapack

mex -R2018a CXXFLAGS='$CXXFLAGS -Wall' ...
    CXXLIBS='$CXXLIBS' ...
    CXXOPTIMFLAGS='$CXXOPTIMFLAGS -O3 -march=native' ...
    reconstruction_LLS2.cpp -lmwlapack

mex -R2018a CXXFLAGS='$CXXFLAGS -Wall' ...
    CXXLIBS='$CXXLIBS' ...
    CXXOPTIMFLAGS='$CXXOPTIMFLAGS -O3 -march=native' ...
    reconstruction_LLS3.cpp -lmwlapack
