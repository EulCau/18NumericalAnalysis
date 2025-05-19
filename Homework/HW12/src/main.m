N = [10, 20, 40, 80, 160];

a1 = 0; b1 = pi / 2;
alpha1 = 3; beta1 = 7;
u1 = @(t) 0; v1 = @(t) -1; w1 = @(t) 0;

a2 = 0; b2 = 1;
alpha2 = 2; beta2 = exp(1) + cos(1);
u2 = @(t) 2 * exp(t); v2 = @(t) -1; w2 = @(t) 0;

errors1 = zeros(size(N)); orders1 = NaN(size(N));
errors2 = zeros(size(N)); orders2 = NaN(size(N));

fprintf("=== BVP 1 ===\n");
fprintf("%6s %12s %10s\n", 'N', 'Error', 'Order');

for i =1:length(N)
    n = N(i);
    ts = linspace(a1, b1, n+1);

    xs = solve_bvp_fd(a1, b1, alpha1, beta1, n, u1, v1, w1, n == N(end));
    xReals = 7 * sin(ts) + 3 * cos(ts);
    errors1(i) = max(abs(xReals - xs));

    if i ~= 1
        orders1(i) = log(errors1(i-1)/errors1(i)) / log(N(i)/N(i-1));
    end

    fprintf("%6d %12.4e %10.4f\n", N(i), errors1(i), orders1(i));
end

fprintf("\n=== BVP 2 ===\n");
fprintf("%6s %12s %10s\n", 'N', 'Error', 'Order');

for i = 1:length(N)
    n = N(i);
    ts = linspace(a2, b2, n+1);

    xs = solve_bvp_fd(a2, b2, alpha2, beta2, n, u2, v2, w2, n == N(end));
    xReals = exp(ts) + cos(ts);
    errors2(i) = max(abs(xReals - xs));

    if i ~= 1
        orders2(i) = log(errors2(i-1)/errors2(i)) / log(N(i)/N(i-1));
    end

    fprintf("%6d %12.4e %10.4f\n", N(i), errors2(i), orders2(i));
end
