function [B, A] = peakFilter(fc, Q, fs, gain)
    A = 10^( gain / 40);
    omega = 2 * pi * fc / fs;
    sn = sin(omega);
    cs = cos(omega);
    alpha = sn / (2.0 * Q);
    %beta = sqrt(A + A);

    b0 = 1 + (alpha * A);
    b1 = -2 * cs;
    b2 = 1 - (alpha * A);
    a0 = 1 + (alpha / A);
    a1 = -2 * cs;
    a2 = 1 - (alpha / A);

    b0 = b0 / a0;
    b1 = b1 / a0;
    b2 = b2 / a0;
    a1 = a1 / a0;
    a2 = a2 / a0;
    a0 = a0 / a0;
    %SOS = [b0,b1,b2,a1,a2];
    B = [b0,b1,b2];
    A = [a0,a1,a2];
    %fvtool([b0,b1,b2],[a0,a1,a2]) %Visualize the filter
end