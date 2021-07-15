#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <vector>
#include <algorithm>
#include "mex.h"
#include "matrix.h"
#include "polysoup.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
// function [mask, p, mnv] = polysoup_merge_vertices_cleanup(polysoup);
{
	// argomenti in ingresso
    Polysoup polysoup (prhs[0]);

    // argomenti in uscita
    plhs[0] = mxCreateLogicalMatrix(polysoup.np,1);
    bool *mask = mxGetLogicals(plhs[0]);
    if (mask == nullptr) mexErrMsgIdAndTxt("MEX:nullptr", "mask is NULL");

    plhs[1] = mxDuplicateArray(mxGetField(prhs[0],0,"p"));
    uint32_t *p = mxGetUint32s(plhs[1]);
    if (p == nullptr) mexErrMsgIdAndTxt("MEX:nullptr", "p is NULL");

    // rimuovi i poligoni con meno di tre vertici distinti
    uint32_t mnv = 3;
    std::vector<uint32_t> unique_vertices;
    for (uint32_t i = 1; i <= polysoup.np; i++) {
        unique_vertices.clear();
        for (uint32_t j = 0; j < polysoup.mnv; j++) {
            uint32_t v = p[i-1+j*polysoup.np];
            if (v == 0) break;
            auto search = std::find(unique_vertices.begin(), unique_vertices.end(), v);
            if (search == std::end(unique_vertices)) {
                unique_vertices.push_back(v);
            }
        }
        uint32_t l = unique_vertices.size();
        if (l >= 3) {
            mask[i-1] = true;
            for (uint32_t j = 0; j < polysoup.mnv; j++) {
                if (j < l) {
                    p[i-1+j*polysoup.np] = unique_vertices[j];
                } else {
                    p[i-1+j*polysoup.np] = 0;
                }
            }
            mnv = std::max(mnv, l);
        }
    }
    plhs[2] = mxCreateDoubleScalar((double)mnv);
}
