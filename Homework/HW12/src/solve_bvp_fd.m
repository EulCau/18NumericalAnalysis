function x = solve_bvp_fd(a, b, alpha, beta, n, u, v, w, show)
    % solve_bvp_fd: 使用有限差分法求解两点边值问题
    % 输入:
    %   a, b: 区间端点
    %   alpha, beta: 边界条件, x(a) = alpha, x(b) = beta
    %   n: 网格数量 (会产生 n + 1 个节点)
    %   u, v, w: 函数句柄, 例如 @(t) sin(t), x'' = u(t) + v(t) x + w(t) x'
    %   show: 控制是否绘制函数图像, 默认会进行绘图
    %
    % 输出:
    %   x: 含边值在内的解向量 x(1) = x(a), x(end) = x(b)
    if nargin < 9
        show = true;  % 默认开启绘图
    end

    h = (b - a) / n;              % 步长
    t = linspace(a, b, n+1);     % 网格节点 t_0, ..., t_n

    % 初始化稀疏三对角系数矩阵 A 和右端项向量 F
    A = zeros(n-1, n-1);  % 系数矩阵, 内部点数量为 n - 1
    F = zeros(n-1, 1);    % 右端向量

    for i = 1:n-1
        ti = t(i+1); % t_1 到 t_{n-1}

        % 近似 x'' ≈ (x_{i-1} - 2x_i + x_{i+1}) / h^2
        % 近似 x' ≈ (x_{i+1} - x_{i-1}) / (2h)

        % 系数
        a_i = 1/h^2 - w(ti)/(2*h);
        b_i = -2/h^2 - v(ti);
        c_i = 1/h^2 + w(ti)/(2*h);

        % 设置三对角矩阵中的值
        if i > 1
            A(i, i-1) = a_i;
        else
            F(i) = F(i) - a_i * alpha;
        end

        A(i, i) = b_i;

        if i < n-1
            A(i, i+1) = c_i;
        else
            F(i) = F(i) - c_i * beta;
        end

        % 添加 u(ti)
        F(i) = F(i) + u(ti);
    end

    % 求解线性方程组
    x_inner = A \ F;

    % 构造包含边值的完整解
    x = [alpha; x_inner; beta]';

    % 可视化结果
    if  show
        figure;
        plot(t, x, '-');
        xlabel('t'); ylabel('x(t)');
        title('解的有限差分法近似');
        grid on;
    end
end
