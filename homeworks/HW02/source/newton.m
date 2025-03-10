% 定义函数 f(x)
f = @(x) 1 ./ (1 + 25 * x.^2);

% 定义误差计算函数
compute_error = @(f, p, y) max(abs(f(y) - p(y)));

% 定义插值节点生成函数
generate_nodes = @(N, type) ...
    strcmp(type, 'equidistant') * (1 - 2 * (0:N) / N) + ...
    strcmp(type, 'chebyshev') * (-cos((2 * (0:N) + 1) * pi / (2 * N + 2)));

% 计算并比较误差
N_values = [5, 10, 20, 40];
y = (0:100) / 50 - 1; % y_i = i/50 - 1, i = 0, 1, ..., 100

for k = 1:length(N_values)
    N = N_values(k);
    
    % 等距节点
    x_equidistant = generate_nodes(N, 'equidistant');
    y_equidistant = f(x_equidistant);
    p_equidistant = @(xi) arrayfun(@(x) newton_interpolation(x_equidistant, y_equidistant, x), xi);

    errorsEqu = compute_error(f, p_equidistant, y);
    
    % 切比雪夫节点
    x_chebyshev = generate_nodes(N, 'chebyshev');
    y_chebyshev = f(x_chebyshev);
    p_chebyshev = @(xi) arrayfun(@(x) newton_interpolation(x_chebyshev, y_chebyshev, x), xi);
    errorsChe = compute_error(f, p_chebyshev, y);
    
    % 输出误差比较结果
    fprintf('N=%d\n', N);
    fprintf('Max Error of grid (1) : %.10f\n', errorsEqu);
    fprintf('Max Error of grid (2) : %.10f\n', errorsChe);

    % 绘制 N = 10 时的函数图像
    if k == 3
        x_plot = linspace(-0.8, 0.8, 100);
        y_plot = f(x_plot);
        figure;
        plot(x_plot, y_plot, 'k', 'LineWidth', 1.5); hold on;
        plot(x_plot, p_equidistant(x_plot), 'r--', 'LineWidth', 1);
        plot(x_plot, p_chebyshev(x_plot), 'b-.', 'LineWidth', 1);
        legend('f(x)', '等距节点插值', '切比雪夫节点插值');
        title('N = 20 时的函数插值比较, -0.8 <= x <= 0.8');
        xlabel('x');
        ylabel('y');
        grid on;
        
        x_plot = linspace(0.8, 1, 100);
        y_plot = f(x_plot);
        figure;
        plot(x_plot, y_plot, 'k', 'LineWidth', 1.5); hold on;
        plot(x_plot, p_equidistant(x_plot), 'r--', 'LineWidth', 1);
        plot(x_plot, p_chebyshev(x_plot), 'b-.', 'LineWidth', 1);
        legend('f(x)', '等距节点插值', '切比雪夫节点插值');
        title('N = 20 时的函数插值比较, 0.8 <= x <= 1');
        xlabel('x');
        ylabel('y');
        grid on;
    end
end

% 定义计算牛顿插值多项式的函数
function p = newton_interpolation(x_nodes, y_nodes, x_eval)
    N = length(x_nodes);
    divided_diff = y_nodes;
    for j = 2:N
        for i = N:-1:j
            divided_diff(i) = (divided_diff(i) - divided_diff(i-1)) / (x_nodes(i) - x_nodes(i-j+1));
        end
    end
    
    p = divided_diff(N) * ones(size(x_eval));
    for k = N-1:-1:1
        p = divided_diff(k) + (x_eval - x_nodes(k)) .* p;
    end
end
