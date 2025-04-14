function test(f, x, h, M)
    D = richardson_derivative(f, x, h, M);
    fprintf('For f(x) = %s, x = %.2f, M = %d:\n', func2str(f), x, M);
    disp(D);
end

f = @(x) log(x);
x = 3;
h = 1;
M = 3;
test(f, x, h, M);

f = @(x) tan(x);
x = asin(0.8);
h = 1;
M = 4;
test(f, x, h, M);

f = @(x) sin(x^2 + (1/3)*x);
x = 0;
h = 1;
M = 5;
test(f, x, h, M);
