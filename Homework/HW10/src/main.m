f = @(t,x) (t - exp(-t)) / (x + exp(x));
f_exact = @(x,t) x^2 - t^2 + 2 * exp(x) - 2 * exp(-t);
t0 = 0.5;
x0 = fzero(@(x) f_exact(t0, x), -0.63);
t_ref = 1;
x_ref = fzero(@(x) f_exact(t_ref, x), -0.16);
adams_bashforth5_test(f, t0, x0, t_ref, x_ref);
