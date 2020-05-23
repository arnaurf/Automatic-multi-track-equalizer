#include <stdio.h>
#include <vector>
#include "BiquadFilter.h"


int offset(int x, int y, int z, int xSize, int ySize);

float EMA(float x, float y0, float alpha);

void hanning(std::vector<double>& buffer, int size);

std::vector<double> filter2(double b[], double a[], const float* X, int nCoefs, int sizeX);

void filter3(float* output, BiquadFilter filter, const float* X, int sizeX);

float rms(std::vector<double> x, int n);

float* max(std::vector<float> x, int size);

float clamp(float value, float min, float max);

