% 主程序
% 定义函数 f(x) = sin(x)
f = @(x) sin(x);

function test(f, a, b, exact)
    % 设定 N 的不同值
    k_max = 12;  % k的最大值
    N_values = 2.^(1:k_max);  % N = 2^k

    % 存储误差
    error_trap = zeros(1, k_max);
    error_simpson = zeros(1, k_max);

    fprintf('a = %.4f, b = %.4f\n', a, b);

    for k = 1:k_max
        N = 2^k;
    
        % 计算复化梯形法的积分
        I_trap = trapezoidal_rule(f, a, b, N);
    
        % 计算复化Simpson法的积分
        I_simpson = simpson_rule(f, a, b, N);
    
        % 计算误差
        error_trap(k) = abs(I_trap - exact);
        error_simpson(k) = abs(I_simpson - exact);
    
        % 打印结果
        fprintf('    k = %d, N = %d\n', k, N);
        fprintf('    Trapezoidal: I = %.4f  Error = %.4e\n', I_trap, error_trap(k));
        fprintf('    Simpson:     I = %.4f  Error = %.4e\n\n', I_simpson, error_simpson(k));
    end

    % 计算收敛阶
    error_trap_old = error_trap(1:end-1);
    error_trap_new = error_trap(2:end);
    N_trap_old = N_values(1:end-1);
    N_trap_new = N_values(2:end);

    Ord_trap = -log(error_trap_old ./ error_trap_new) ./ log(N_trap_old ./ N_trap_new);

    error_simpson_old = error_simpson(1:end-1);
    error_simpson_new = error_simpson(2:end);
    N_simpson_old = N_values(1:end-1);
    N_simpson_new = N_values(2:end);

    Ord_simpson = -log(error_simpson_old ./ error_simpson_new) ./ log(N_simpson_old ./ N_simpson_new);

    % 输出收敛阶
    fprintf('Order (Trapezoidal):\n');
    disp(Ord_trap);
    fprintf('Order (Simpson):\n');
    disp(Ord_simpson);
end

% 设置积分区间
a1 = 0; b1 = 4;   % 第一个积分区间
exact1 = -cos(b1) + cos(a1); % 积分结果的准确值

a2 = 0; b2 = 2*pi; % 第二个积分区间
exact2 = 0; % 积分结果的准确值

test(f, a1, b1, exact1)
test(f, a2, b2, exact2)
