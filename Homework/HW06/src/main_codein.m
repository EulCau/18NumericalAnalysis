f = @(x) 1 ./ (1 + 25 * x.^2);
exact = 0.4 * atan(5);

fprintf('等距节点结果：\n');
fprintf('N\tPolyIntegral\tFuncIntegral\tError\n');

for N = 5:5:40
    x_nodes = linspace(-1, 1, N + 1);
    polyInt = lagrange_integral(f, x_nodes, -1, 1);
    err = abs(polyInt - exact);
    fprintf('%d\t%.8f\t%.8f\t%.2e\n', N, polyInt, exact, err);
end

fprintf('\nChebyshev 节点结果：\n');
fprintf('N\tPolyIntegral\tFuncIntegral\tError\n');

for N = 5:5:40
    x_nodes = -cos((1:N+1) * pi / (N + 2));
    polyInt = lagrange_integral(f, x_nodes, -1, 1);
    err = abs(polyInt - exact);
    fprintf('%d\t%.8f\t%.8f\t%.2e\n', N, polyInt, exact, err);
end
