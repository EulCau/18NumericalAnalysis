N_list = [5, 10, 20, 40];
error_linear = zeros(size(N_list));
error_cubic = zeros(size(N_list));

for k = 1:length(N_list)
    N = N_list(k);
    h = 1 / N;
    x = linspace(0, 1, N + 1);
    f = exp(x);
    midpoints = x(1:end-1) + h / 2;
    
    % 分片线性样条的误差计算
    s_linear = (f(1:end-1) + f(2:end)) / 2;
    error_linear(k) = max(abs(exp(midpoints) - s_linear));
    
    % 三次样条构造
    A = zeros(N + 1, N + 1);
    b = zeros(N + 1, 1);
    
    % 左边界条件
    A(1, 1) = h / 3;
    A(1, 2) = h / 6;
    b(1) = (f(2) - f(1)) / h - 1;
    
    % 右边界条件
    A(end, end-1) = h / 6;
    A(end, end) = h / 3;
    b(end) = exp(1) - (f(end) - f(end-1)) / h;
    
    % 内部节点方程
    for i = 2:N
        j = i - 1;
        A(i, j) = 1;
        A(i, j + 1) = 4;
        A(i, j + 2) = 1;
        b(i) = (6 / h^2) * (f(j + 2) - 2 * f(j + 1) + f(j));
    end
    
    M = A \ b; % 解方程
    
    % 计算三次样条在中点的误差
    s_cubic = (f(1:N) + f(2:N+1)) / 2 - ((M(1:N) + M(2:N+1))') * h^2 / 48;
    error_cubic(k) = max(abs(exp(midpoints) - s_cubic));
end

% 计算收敛阶
ord_linear = zeros(1, length(N_list));
ord_cubic = zeros(1, length(N_list));

for k = 2:length(N_list)
    ord_linear(k) = log(error_linear(k-1)/error_linear(k)) / log(N_list(k)/N_list(k-1));
    ord_cubic(k) = log(error_cubic(k-1)/error_cubic(k)) / log(N_list(k)/N_list(k-1));
end

% 显示结果
fprintf('N\t线性样条误差\t收敛阶\t三次样条误差\t收敛阶\n');
for k = 1:length(N_list)
    fprintf('%d\t%.4e\t%.4f\t%.4e\t%.4f\n', N_list(k), error_linear(k), ord_linear(k), error_cubic(k), ord_cubic(k));
end