f = @(t,x) (t - exp(-t)) / (x + exp(x));
f_exact = @(t,x) x^2 - t^2 + 2 * exp(x) - 2 * exp(-t) + 1 - 2 * exp(1);
t0 = 0;
x0 = 1;
t_ref = 1;
x_ref = fzero(@(x) f_exact(t_ref, x), 1);
adams_bashforth5(f, t0, x0, t_ref, x_ref);
