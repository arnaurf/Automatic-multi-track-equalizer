#include <stdio.h>
#include <vector>
#include "BiquadFilter.h"

#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif
int offset(int x, int y, int z, int xSize, int ySize) {
    return (z * xSize * ySize) + (y * xSize) + x;
}

float EMA(float x, float y0, float alpha) {
    return (1 - alpha) * x + alpha * y0;
}

void hanning(std::vector<double>& buffer, int size) {
    for (int i = 0; i < size; i++) {
        buffer[i] = 0.5 * (1 - cos(2 * M_PI * i / size));
    }
}

std::vector<double> filter2(double b[], double a[], const float* X, int nCoefs, int sizeX) {

    //float* z = (float*)calloc(sizeX, sizeof(float));
    for (int i = 0; i < nCoefs; i++) {
        b[i] = b[i] / a[0];
        a[i] = a[i] / a[0];
    }

    std::vector<double> Y(sizeX);//(double*)calloc(sizeX, sizeof(double));

    for (int n = 0; n < sizeX; n++) {
        double auxX = 0;
        double auxY = 0;
        for (int m = 0; m < nCoefs; m++) {
            if (n - m >= 0)
                auxX = auxX + X[n - m] * b[m];
        }
        for (int m = 1; m < nCoefs; m++) {
            if (n - m >= 0)
                auxY = auxY + Y[n - m] * a[m];
        }
        Y[n] = (auxX - auxY) / double(a[0]);
    }

    /*for (int m = 1; m < sizeX; m++) {
        Y[m] = b[0] * X[m] + z[0];
        for (int i = 2; i < sizeF; i++) {
            z[i - 1] = b[i] * X[m] + z[i] - a[i] * Y[m];
        }
    }
    //z = z[1:n - 1];*/
    return Y;
}

void filter3(float* output, BiquadFilter filter, const float* X, int sizeX) {
    double* coefs = filter.getCoefs();
    double a[3] = { coefs[0], coefs[1], coefs[2] };
    double b[3] = { coefs[3], coefs[4], coefs[5] };
    int sizeB = 3;
    int sizeA = 3;

    //float* z = (float*)calloc(sizeX, sizeof(float));
    for (int i = 0; i < sizeB; i++) {
        b[i] = b[i] / a[0];
    }
    for (int i = 0; i < sizeA; i++) {
        a[i] = a[i] / a[0];
    }
    std::vector<double> Y(sizeX);//(double*)calloc(sizeX, sizeof(double));

    for (int n = 0; n < sizeX; n++) {
        double auxX = 0;
        double auxY = 0;
        for (int m = 0; m < sizeB; m++) {
            if (n - m >= 0)
                auxX = auxX + X[n - m] * b[m];
        }
        for (int m = 1; m < sizeA; m++) {
            if (n - m >= 0)
                auxY = auxY + Y[n - m] * a[m];
        }
        Y[n] = (auxX - auxY) / double(a[0]);
    }

    for (int i = 0; i < sizeX; i++)
        output[i] = Y[i];

}

float rms(std::vector<double> x, int n)
{
    double sum = 0;

    for (int i = 0; i < n; i++)
        sum += pow(x[i], 2);

    return sqrt(sum / n);
}


//It returns the [value, index] of the maximum value inside of x.
float* max(std::vector<float> x, int size) {

    float temp_max[2] = { x[0],0 };
    for (int i = 0; i < size; i++) {
        if (x[i] > temp_max[0]) {
            temp_max[0] = x[i];
            temp_max[1] = i;
        }
    }
    return temp_max;
}

float clamp(float value, float min, float max) {
    if (value > max)
        return max;
    if (value < min)
        return min;
    return value;
}