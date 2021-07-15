#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <vector>
#include <unordered_map>
#include "mex.h"
#include "matrix.h"
#include "polysoup.h"

struct Edges {
	uint32_t ne;
	std::vector<uint32_t> type;
    std::vector<uint32_t> v1;
    std::vector<uint32_t> v2;
    std::vector<uint32_t> cp;
    std::vector<uint32_t> cm;
};

struct Cells {
    uint32_t nc;
    uint32_t mne;
    uint32_t *ne;
    int32_t *e;
};

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
// function [edges, cells] = polymesh_from_polysoup_merge_edges(polysoup);
{
	// argomenti in ingresso
    Polysoup polysoup (prhs[0]);

    // argomenti in uscita
    const char *fieldnames_edges[6] = {"ne", "type", "v1", "v2", "cp", "cm"};
    const char *fieldnames_cells[4] = {"nc", "mne", "ne", "e"};
    plhs[0] = mxCreateStructMatrix(1, 1, 6, fieldnames_edges);
    plhs[1] = mxCreateStructMatrix(1, 1, 4, fieldnames_cells);

    // inizializza lo struct "cells"
    Cells cells;

    mxArray *cells_nc = mxCreateDoubleScalar(polysoup.np);
    mxSetField(plhs[1], 0, "nc", cells_nc);
    cells.nc = (uint32_t)mxGetScalar(cells_nc);

    mxArray *cells_mne = mxCreateDoubleScalar(polysoup.mnv);
    mxSetField(plhs[1], 0, "mne", cells_mne);
    cells.mne = (uint32_t)mxGetScalar(cells_mne);

    mxArray *cells_ne = mxCreateNumericMatrix(cells.nc,1,mxUINT32_CLASS,mxREAL);
    mxSetField(plhs[1], 0, "ne", cells_ne);
    cells.ne = mxGetUint32s(cells_ne);
    if (cells.ne == nullptr) mexErrMsgIdAndTxt("MEX:nullptr", "cells.ne is NULL");

    mxArray *cells_e = mxCreateNumericMatrix(cells.nc,cells.mne,mxINT32_CLASS,mxREAL);
    mxSetField(plhs[1], 0, "e", cells_e);
    cells.e = mxGetInt32s(cells_e);
    if (cells.e == nullptr) mexErrMsgIdAndTxt("MEX:nullptr", "cells.e is NULL");

    // ricostruisci gli spigoli
    Edges edges;
    std::unordered_map<uint64_t,uint32_t> hash_table;
    uint32_t counter_e = 1;
    for (uint32_t i = 1; i <= polysoup.np; i++) {
        uint32_t ne = 0;
        for (uint32_t j = 0; j < polysoup.mnv; j++) {
            if (polysoup.p[i-1+j*polysoup.np] != 0) {
                ne += 1;
            }
        }
        cells.ne[i-1] = ne;
        for (uint32_t j = 0; j < ne; j++) {
            uint32_t jnext = (j+1) % ne;
            uint32_t v1 = polysoup.p[i-1+j*polysoup.np];
            uint32_t v2 = polysoup.p[i-1+jnext*polysoup.np];
            uint64_t v12;
            if (v1 < v2) {
                v12 = ((uint64_t)v1 << 32) | (uint64_t)v2;
            } else {
                v12 = ((uint64_t)v2 << 32) | (uint64_t)v1;
            }
            auto search = hash_table.find(v12);
            if (search == hash_table.end()) {
                // non trovato
                hash_table.insert({v12,counter_e});
                edges.type.push_back(1);
                edges.v1.push_back(v1);
                edges.v2.push_back(v2);
                edges.cp.push_back(i);
                edges.cm.push_back(0);
                cells.e[i-1+j*cells.nc] = counter_e;
                counter_e += 1;
            } else {
                // trovato
                uint32_t k = search->second;
                edges.type[k-1] += 1;
                edges.cm[k-1] = i;
                cells.e[i-1+j*cells.nc] = -(int32_t)k;
            }
        }
    }
    edges.ne = counter_e - 1;

    // copia i vettori dello struct Edges nei campi di plhs[0]
    mxArray *edges_ne = mxCreateDoubleScalar((double)edges.ne);
    mxSetField(plhs[0], 0, "ne", edges_ne);

    mxArray *edges_type = mxCreateNumericMatrix(edges.ne,1,mxUINT32_CLASS,mxREAL);
    mxSetField(plhs[0], 0, "type", edges_type);
    memcpy(mxGetUint32s(edges_type), edges.type.data(), 4*edges.ne);

    mxArray *edges_v1 = mxCreateNumericMatrix(edges.ne,1,mxUINT32_CLASS,mxREAL);
    mxSetField(plhs[0], 0, "v1", edges_v1);
    memcpy(mxGetUint32s(edges_v1), edges.v1.data(), 4*edges.ne);

    mxArray *edges_v2 = mxCreateNumericMatrix(edges.ne,1,mxUINT32_CLASS,mxREAL);
    mxSetField(plhs[0], 0, "v2", edges_v2);
    memcpy(mxGetUint32s(edges_v2), edges.v2.data(), 4*edges.ne);

    mxArray *edges_cp = mxCreateNumericMatrix(edges.ne,1,mxUINT32_CLASS,mxREAL);
    mxSetField(plhs[0], 0, "cp", edges_cp);
    memcpy(mxGetUint32s(edges_cp), edges.cp.data(), 4*edges.ne);

    mxArray *edges_cm = mxCreateNumericMatrix(edges.ne,1,mxUINT32_CLASS,mxREAL);
    mxSetField(plhs[0], 0, "cm", edges_cm);
    memcpy(mxGetUint32s(edges_cm), edges.cm.data(), 4*edges.ne);
}
