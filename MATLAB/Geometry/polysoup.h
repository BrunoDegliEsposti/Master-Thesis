#pragma once
#include <cstdint>
#include "mex.h"
#include "matrix.h"

struct Polysoup {
	uint32_t nv;
	double *vx;
	double *vy;
	uint32_t np;
	uint32_t mnv;
	uint32_t *p;
	double *cx;
	double *cy;
	Polysoup(const mxArray *in)
	{
		nv = (uint32_t)mxGetScalar(mxGetField(in,0,"nv"));
		vx = mxGetDoubles(mxGetField(in,0,"vx"));
		vy = mxGetDoubles(mxGetField(in,0,"vy"));
		np = (uint32_t)mxGetScalar(mxGetField(in,0,"np"));
		mnv = (uint32_t)mxGetScalar(mxGetField(in,0,"mnv"));
		p = mxGetUint32s(mxGetField(in,0,"p"));
		cx = mxGetDoubles(mxGetField(in,0,"cx"));
		cy = mxGetDoubles(mxGetField(in,0,"cy"));
		if (vx == nullptr) mexErrMsgIdAndTxt("MEX:nullptr", "polysoup.vx = NULL");
		if (vy == nullptr) mexErrMsgIdAndTxt("MEX:nullptr", "polysoup.vx = NULL");
		if (p == nullptr) mexErrMsgIdAndTxt("MEX:nullptr", "polysoup.vx = NULL");
		if (cx == nullptr) mexErrMsgIdAndTxt("MEX:nullptr", "polysoup.vx = NULL");
		if (cy == nullptr) mexErrMsgIdAndTxt("MEX:nullptr", "polysoup.vx = NULL");
	}
};
