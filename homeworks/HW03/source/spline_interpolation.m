clc; clear; close all;

% 计算误差和收敛阶
N_values = [5, 10, 20, 40];
errors_linear = zeros(length(N_values),1);
errors_cubic = zeros(length(N_values),1);

for idx = 1:length(N_values)
    N = N_values(idx);
    h = 1/N;
    x = linspace(0, 1, N+1);
    y = exp(x);

    % 计算线性样条误差
    x_half = x(1:end-1) + h/2;
    y_half_true = exp(x_half);
    y_half_linear = interp_linear(x, y, x_half);
    errors_linear(idx) = max(abs(y_half_true - y_half_linear));

    % 计算三次样条误差
    y_half_cubic = interp_cubic_spline(x, y, x_half);
    errors_cubic(idx) = max(abs(y_half_true - y_half_cubic));
end

% 计算收敛阶
order_linear = log(errors_linear(1:end-1) ./ errors_linear(2:end)) ./ log(N_values(2:end) ./ N_values(1:end-1))';
order_cubic = log(errors_cubic(1:end-1) ./ errors_cubic(2:end)) ./ log(N_values(2:end) ./ N_values(1:end-1))';

% 输出结果
fprintf('N\t线性样条误差\t收敛阶\t三次样条误差\t收敛阶\n');
fprintf('%d\t%.4e\t%.4f\t%.4e\t%.4f\n', [N_values', errors_linear, [0; order_linear], errors_cubic, [0; order_cubic]]');

% 线性插值函数
function y_interp = interp_linear(x, y, x_half)
    N = length(x) - 1;
    y_interp = zeros(size(x_half));
    for i = 1:N
        idx = (x_half >= x(i)) & (x_half <= x(i+1));
        y_interp(idx) = y(i) + (y(i+1) - y(i)) / (x(i+1) - x(i)) * (x_half(idx) - x(i));
    end
end

% 三次样条插值函数
function y_interp = interp_cubic_spline(x, y, x_half)
    N = length(x) - 1;
    h = diff(x);
    
    % 计算二阶导数矩阵
    A = zeros(N+1, N+1);
    rhs = zeros(N+1, 1);
    
    % 边界条件 S'(0) = 1, S'(1) = exp(1)
    A(1,1) = 1; rhs(1) = 1;
    A(N+1,N+1) = 1; rhs(N+1) = exp(1);

    % 组装三对角矩阵
    for i = 2:N
        A(i, i-1) = h(i-1);
        A(i, i) = 2 * (h(i-1) + h(i));
        A(i, i+1) = h(i);
        rhs(i) = 3 * ((y(i+1) - y(i)) / h(i) - (y(i) - y(i-1)) / h(i-1));
    end

    % 解方程求 c
    c = A \ rhs;

    % 计算 b 和 d
    b = zeros(N,1);
    d = zeros(N,1);
    for i = 1:N
        b(i) = (y(i+1) - y(i)) / h(i) - h(i) * (c(i+1) + 2*c(i)) / 3;
        d(i) = (c(i+1) - c(i)) / (3*h(i));
    end

    % 插值计算
    y_interp = zeros(size(x_half));
    for i = 1:N
        idx = (x_half >= x(i)) & (x_half <= x(i+1));
        dx = x_half(idx) - x(i);
        y_interp(idx) = y(i) + b(i) * dx + c(i) * dx.^2 + d(i) * dx.^3;
    end
end
