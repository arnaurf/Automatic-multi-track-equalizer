
#include <math.h>
#include "BiquadFilter.h"

BiquadFilter::BiquadFilter(float fc, float Q, float fs, float gain) {

	double omega = 2 * M_PI * fc / fs;
	double sn = sin(omega);
	double cs = cos(omega);
	double alpha = sn / double(2.0f * Q);
	double  A = pow(10, gain / 40);
	double beta = sqrt(A + A);
	b0 = 1 + (alpha * A);
	b1 = -2 * cs;
	b2 = 1 - (alpha * A);
	a0 = 1 + (alpha / A);
	a1 = -2 * cs;
	a2 = 1 - (alpha / A);

	b0 = b0 / a0;
	b1 = b1 / a0;
	b2 = b2 / a0;
	a0 = a0 / a0;
	a1 = a1 / a0;
	a2 = a2 / a0;

}

double* BiquadFilter::getCoefs() {
	double aux[6] = { a0,a1,a2,b0,b1,b2 };
	return aux;
}