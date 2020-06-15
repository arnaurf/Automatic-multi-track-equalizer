function b=bark(f)
% b=bark(f)
% Converts frequency to bark scale
% Frequency should be specified in Hertz

b = 13*atan(0.76*f/1000) + 3.5*atan((f/7500).^2); 

return;