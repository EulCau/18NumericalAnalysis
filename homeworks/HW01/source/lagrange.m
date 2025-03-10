% 定义函数 f(x)
f = @(x) 1 ./ (1 + x.^2);

% 定义误差计算函数
compute_error = @(f, p, y) max(abs(f(y) - p(y)));

% 定义插值节点生成函数
generate_nodes = @(N, type) ...
    strcmp(type, 'equidistant') * (5 - 10 * (0:N) / N) + ...
    strcmp(type, 'chebyshev') * (5 * cos((2 * (0:N) + 1) * pi / (2 * N + 2)));

% 定义拉格朗日插值函数
lagrange_interpolation = @(x, y, xi) ...
    sum(arrayfun(@(i) y(i+1) * ...
    prod((xi - x([1:i, i+2:end])) ./ ...
    (x(i+1) - x([1:i, i+2:end]))), 0:length(x)-1));

% 计算并比较误差
N_values = [5, 10, 20, 40];
errors = zeros(length(N_values), 2);
y = (-50:50) / 10; % y_i = i/10 - 5, i = 0, 1, ..., 100

for k = 1:length(N_values)
    N = N_values(k);
    
    % 等距节点
    x_equidistant = generate_nodes(N, 'equidistant');
    y_equidistant = f(x_equidistant);
    p_equidistant = @(xi) arrayfun(@(x) lagrange_interpolation(x_equidistant, y_equidistant, x), xi);
    errors(k, 1) = compute_error(f, p_equidistant, y);
    
    % 切比雪夫节点
    x_chebyshev = generate_nodes(N, 'chebyshev');
    y_chebyshev = f(x_chebyshev);
    p_chebyshev = @(xi) arrayfun(@(x) lagrange_interpolation(x_chebyshev, y_chebyshev, x), xi);
    errors(k, 2) = compute_error(f, p_chebyshev, y);
    
    % 输出误差比较结果
    fprintf('N=%d\n', N);
    fprintf('Max Error of grid (1) : %.10f\n', errors(k, 1));
    fprintf('Max Error of grid (2) : %.10f\n', errors(k, 2));
    
    % 绘制 N = 10 时的函数图像
    if k == 2
        N_plot = 10;
        x_plot = linspace(-5, 5, 1000);
        y_plot = f(x_plot);
        figure;
        plot(x_plot, y_plot, 'k', 'LineWidth', 1.5); hold on;
        plot(x_plot, p_equidistant(x_plot), 'r--', 'LineWidth', 1);
        plot(x_plot, p_chebyshev(x_plot), 'b-.', 'LineWidth', 1);
        legend('f(x)', '等距节点插值', '切比雪夫节点插值');
        title('N = 10 时的函数插值比较');
        xlabel('x');
        ylabel('y');
        grid on;
    end
end
