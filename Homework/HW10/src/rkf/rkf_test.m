f = @(x, y) exp(y*x) - cos(y - x);
x0 = 1;
y0 = 3;
rkf_adaptive_solver(f, x0, y0, true);
