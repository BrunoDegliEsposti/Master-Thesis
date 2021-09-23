#pragma once
#include <cstdint>
#include "mex.h"
#include "matrix.h"

struct Method {
	uint32_t order;
	char least_squares_type;
	double WENO_epsilon;
	double WENO_power;
	Method(const mxArray *in)
	{
		mxArray *order_field = mxGetField(in,0,"order");
		if (order_field == nullptr) mexErrMsgIdAndTxt("MEX:nullptr", "method has no field order");
		order = (uint32_t)mxGetScalar(order_field);

		mxArray *least_squares_type_field = mxGetField(in,0,"least_squares_type");
		if (least_squares_type_field != nullptr) {
			char *least_squares_type_array = mxArrayToUTF8String(least_squares_type_field);
			if (least_squares_type_array != nullptr) {
				least_squares_type = least_squares_type_array[0];
			} else {
				mexErrMsgIdAndTxt("MEX:nullptr", "method.least_squares_type is not a char");
			}
		} else {
			least_squares_type = 'P';
		}

		mxArray *WENO_epsilon_field = mxGetField(in,0,"WENO_epsilon");
		if (WENO_epsilon_field != nullptr) {
			WENO_epsilon = mxGetScalar(WENO_epsilon_field);
		} else {
			WENO_epsilon = 1e-5;
		}

		mxArray *WENO_power_field = mxGetField(in,0,"WENO_power");
		if (WENO_power_field != nullptr) {
			WENO_power = mxGetScalar(WENO_power_field);
		} else {
			WENO_power = 4.0;
		}
	}
};
