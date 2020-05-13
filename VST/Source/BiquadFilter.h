#pragma once
#define M_PI 3.14159265358979323846

class BiquadFilter {

public:
	double a0, a1, a2, b0, b1, b2;
	BiquadFilter(float fc, float Q, float fs, float gain);
	double* getCoefs();
};