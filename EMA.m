%x -> new value     y0 -> last value
function y = EMA(x, y0, alpha)

y = (1-alpha)*x + alpha*y0;

end