f = @(t,x) (t - exp(-t)) / (x + exp(x));
x_ref = fzero(@(x) x^2 - 1 + 2*exp(x) - 2*exp(-1), 0.5);
adams_bashforth5_test(f, x_ref);