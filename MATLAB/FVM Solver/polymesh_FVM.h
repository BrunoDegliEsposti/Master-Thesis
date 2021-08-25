#pragma once
#include <cstdint>
#include "mex.h"
#include "matrix.h"
#include "../Geometry/polymesh.h"

struct VerticesFVM: Vertices {
	VerticesFVM(const mxArray *in): Vertices(in) {}
};

struct EdgesFVM: Edges {

    // FVM attributes
    uint32_t nq;
    double *qx;
    double *qw;
    double *up;
    double *um;
    double *tnf;
    double *mws;

	EdgesFVM(const mxArray *in): Edges(in)
	{
		nq = (uint32_t)mxGetScalar(mxGetField(in,0,"nq"));
		qx = mxGetDoubles(mxGetField(in,0,"qx"));
		qw = mxGetDoubles(mxGetField(in,0,"qw"));
		up = mxGetDoubles(mxGetField(in,0,"up"));
		um = mxGetDoubles(mxGetField(in,0,"um"));
		tnf = mxGetDoubles(mxGetField(in,0,"tnf"));
		mws = mxGetDoubles(mxGetField(in,0,"mws"));
		if (qx == nullptr) mexErrMsgIdAndTxt("MEX:nullptr", "edges.qx is NULL");
		if (qw == nullptr) mexErrMsgIdAndTxt("MEX:nullptr", "edges.qw is NULL");
		// if (up == nullptr) mexErrMsgIdAndTxt("MEX:nullptr", "edges.up is NULL");
		// if (um == nullptr) mexErrMsgIdAndTxt("MEX:nullptr", "edges.um is NULL");
		// if (tnf == nullptr) mexErrMsgIdAndTxt("MEX:nullptr", "edges.tnf is NULL");
		// if (mws == nullptr) mexErrMsgIdAndTxt("MEX:nullptr", "edges.mws is NULL");
	}
};

struct CellsFVM: Cells {

    // FVM attributes
    uint32_t nu;
    double *u;
    double *mws;
    double *camb;

    CellsFVM(const mxArray *in): Cells(in)
    {
    	nu = (uint32_t)mxGetScalar(mxGetField(in,0,"nu"));
    	u = mxGetDoubles(mxGetField(in,0,"u"));
    	mws = mxGetDoubles(mxGetField(in,0,"mws"));
    	camb = mxGetDoubles(mxGetField(in,0,"camb"));
    }
};

struct Vec2 {
	double x;
	double y;	
};

Vec2 edge_lerp(double t, Vertices &vertices, Edges &edges, uint32_t j)
{
	uint32_t v1 = edges.v1[j-1];
	uint32_t v2 = edges.v2[j-1];
	double v1x = vertices.x[v1-1];
	double v2x = vertices.x[v2-1];
	double v1y = vertices.y[v1-1];
	double v2y = vertices.y[v2-1];
	return {(1.0-t)*v1x + t*v2x,
		    (1.0-t)*v1y + t*v2y};
}
