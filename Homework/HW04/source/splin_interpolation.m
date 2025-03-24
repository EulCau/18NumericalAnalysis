% MATLAB 代码实现样条插值并计算误差

node_num = [4, 11, 21, 41, 161];
lin_err = zeros(1, length(node_num));
cub_err = zeros(1, length(node_num));
her_err = zeros(1, length(node_num));

for i = 1:length(node_num)
    nodes = linspace(0, 5, node_num(i));
    errors = splin_interpolation_main(nodes, i == 1);
    lin_err(i) = errors.Linear;
    cub_err(i) = errors.Cubic;
    her_err(i) = errors.Hermite;
end

lin_ord = log(lin_err(1:end-1) ./ lin_err(2:end)) ./ log(node_num(2:end) ./ node_num(1:end-1));
cub_ord = log(cub_err(1:end-1) ./ cub_err(2:end)) ./ log(node_num(2:end) ./ node_num(1:end-1));
her_ord = log(her_err(1:end-1) ./ her_err(2:end)) ./ log(node_num(2:end) ./ node_num(1:end-1));

disp('lin_ord');
disp(lin_ord);
disp('cub_ord');
disp(cub_ord);
disp('her_ord');
disp(her_ord);

function errors = splin_interpolation_main(nodes, ifplot)
    values = f(nodes);
    derivatives = -2 * nodes ./ (1 + nodes.^2).^(2);
    grid = linspace(0, 5, 100);
    
    errors.Linear = compute_error(nodes, @linear_spline, values);
    errors.Cubic = compute_error(nodes, @natural_cubic_spline, values);
    errors.Hermite = compute_error(nodes, @hermite_spline, values, derivatives);
    
    if ifplot
        figure;
        hold on;
        plot(grid, f(grid), 'k', 'DisplayName', 'Original Function');
        plot(grid, arrayfun(@(x) linear_spline(nodes, values, x), grid), 'r', 'DisplayName', 'Linear Spline');
        plot(grid, arrayfun(@(x) natural_cubic_spline(nodes, values, x), grid), 'b', 'DisplayName', 'Cubic Spline');
        plot(grid, arrayfun(@(x) hermite_spline(nodes, values, x, derivatives), grid), 'g', 'DisplayName', 'Hermite Spline');
        legend;
        hold off;
    end
    
    fprintf('L1 error when interval_num = %d\n', length(nodes) - 1);
    disp(errors);
end

function y = f(x)
    y = 1 ./ (1 + x.^2);
end

function y = linear_spline(nodes, values, x)
    idx = find(nodes <= x, 1, 'last');
    if idx == length(nodes)
        idx = idx - 1;
    end
    y = values(idx) + (values(idx+1) - values(idx)) / (nodes(idx+1) - nodes(idx)) * (x - nodes(idx));
end

function y = natural_cubic_spline(nodes, values, x)
    n = length(nodes) - 1;
    h = diff(nodes);
    alpha = zeros(1, n);
    for i = 2:n
        alpha(i) = (3./h(i)) .* (values(i+1) - values(i)) - (3./h(i-1)) .* (values(i) - values(i-1));
    end

    l = ones(1, n+1);
    mu = zeros(1, n);
    z = zeros(1, n+1);

    for i = 2:n
        l(i) = 2 * (nodes(i+1) - nodes(i-1)) - h(i-1) .* mu(i-1);
        mu(i) = h(i) ./ l(i);
        z(i) = (alpha(i) - h(i-1) .* z(i-1)) ./ l(i);
    end

    b = zeros(1, n);
    c = zeros(1, n+1);
    d = zeros(1, n);

    for j = n:-1:1
        c(j) = z(j) - mu(j) .* c(j+1);
        b(j) = (values(j+1) - values(j)) ./ h(j) - h(j) .* (c(j+1) + 2.*c(j)) ./ 3;
        d(j) = (c(j+1) - c(j)) ./ (3.*h(j));
    end

    idx = find(nodes <= x, 1, 'last');
    if idx == length(nodes)
        idx = idx - 1;
    end
    dx = x - nodes(idx);
    y = values(idx) + b(idx).*dx + c(idx).*dx.^2 + d(idx).*dx.^3;
end

function y = hermite_spline(nodes, values, xq, derivatives)
    n = length(nodes) - 1;
    y = zeros(size(xq));

    for i = 1:n
        idx = (xq >= nodes(i) & xq <= nodes(i+1));
        h = nodes(i+1) - nodes(i);
        t = (xq(idx) - nodes(i)) ./ h;

        h00 = (1 + 2 .* t) .* (1 - t).^2;
        h10 = t .* (1 - t).^2;
        h01 = t.^2 .* (3 - 2 .* t);
        h11 = t.^2 .* (t - 1);

        y(idx) = h00 .* values(i) + h10 .* h .* derivatives(i) + h01 .* values(i+1) + h11 .* h .* derivatives(i+1);
    end
end

function err = compute_error(nodes, method, values, varargin)
    midpoints = (nodes(1:end-1) + nodes(2:end)) ./ 2;
    err = max(abs(f(midpoints) - arrayfun(@(x) method(nodes, values, x, varargin{:}), midpoints)));
end
