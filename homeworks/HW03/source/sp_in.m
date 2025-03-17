% 定义函数
f = @(x) exp(x);

% 不同节点数量
N_values = [5, 10, 20, 40];

% 存储误差
errors_linear = [];
errors_cubic = [];

for N = N_values
    % 等距节点
    x_nodes = linspace(0, 1, N+1);
    y_nodes = f(x_nodes);

    % 中点
    x_mid = (x_nodes(1:end-1) + x_nodes(2:end)) / 2;
    y_true = f(x_mid);

    % 一次分片线性样条
    y_linear = zeros(1, length(x_mid));
    for i = 1:N
        % 线性插值公式
        slope = (y_nodes(i+1) - y_nodes(i)) / (x_nodes(i+1) - x_nodes(i));
        intercept = y_nodes(i) - slope * x_nodes(i);
        y_linear(i) = slope * x_mid(i) + intercept;
    end
    errors_linear = [errors_linear, max(abs(y_true - y_linear))];

    % 满足边界条件的三次样条插值
    y_cubic = zeros(1, length(x_mid)); % 用于存储中点的三次样条值

    % 构建三次样条的系数矩阵
    h = diff(x_nodes); % 每段的步长
    A = zeros(N+1, N+1); % 系数矩阵
    b = zeros(N+1, 1); % 右侧常数项

    % 填充系数矩阵和右侧项
    for i = 2:N
        A(i, i-1) = h(i-1);
        A(i, i) = 2 * (h(i-1) + h(i));
        A(i, i+1) = h(i);
        b(i) = 3 * ((y_nodes(i+1) - y_nodes(i)) / h(i) - (y_nodes(i) - y_nodes(i-1)) / h(i-1));
    end

    % 边界条件: S'(0) = 1, S'(1) = exp(1)
    A(1,1) = 2*h(1);
    A(1,2) = h(1);
    b(1) = 3 * ((y_nodes(2) - y_nodes(1)) / h(1) - 1);

    A(N+1,N) = h(N);
    A(N+1,N+1) = 2*h(N);
    b(N+1) = 3 * (exp(1) - (y_nodes(N+1) - y_nodes(N)) / h(N));

    % 求解方程组，得到每个节点处的二阶导数 M
    M = A \ b;

    % 计算每段的三次样条系数
    for i = 1:N
        a = y_nodes(i);
        b = (y_nodes(i+1) - y_nodes(i)) / h(i) - h(i) * (2 * M(i) + M(i+1)) / 3;
        c = M(i) / 2;
        d = (M(i+1) - M(i)) / (3 * h(i));
        % 计算中点的样条值
        y_cubic(i) = a + b * (x_mid(i) - x_nodes(i)) + c * (x_mid(i) - x_nodes(i))^2 + d * (x_mid(i) - x_nodes(i))^3;
    end

    % 记录最大误差
    errors_cubic = [errors_cubic, max(abs(y_true - y_cubic))];
end

% 计算收敛阶
orders_linear = - diff(log(errors_linear)) ./ diff(log(N_values));
orders_cubic = - diff(log(errors_cubic)) ./ diff(log(N_values));

disp('一次分片线性样条误差：');
disp(errors_linear);
disp('收敛阶：');
disp(orders_linear);

disp('三次分片线性样条误差：');
disp(errors_cubic);
disp('收敛阶：');
disp(orders_cubic);
