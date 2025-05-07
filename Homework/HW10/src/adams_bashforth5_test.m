function adams_bashforth5_test(f, x_ref)
    % 设置步数
    ks = 3:8;
    Ns = 2.^ks;
    hs = 1 ./ Ns;
    errors = zeros(size(Ns));

    for j = 1:length(Ns)
        N = Ns(j);
        h = hs(j);
        t = 0:h:1;

        % 初始化解向量
        x = zeros(1, N+1);

        % 用 RK4 初始化前 5 步
        for i = 1:4
            k1 = f(t(i), x(i));
            k2 = f(t(i)+h/2, x(i)+h/2*k1);
            k3 = f(t(i)+h/2, x(i)+h/2*k2);
            k4 = f(t(i)+h, x(i)+h*k3);
            x(i+1) = x(i) + (h/6)*(k1 + 2*k2 + 2*k3 + k4);
        end

        % 用 AB5 方法计算后续值
        for i = 5:N
            f1 = f(t(i),   x(i));
            f2 = f(t(i-1), x(i-1));
            f3 = f(t(i-2), x(i-2));
            f4 = f(t(i-3), x(i-3));
            f5 = f(t(i-4), x(i-4));
            x(i+1) = x(i) + (h/720)*(1901*f1 - 2774*f2 + 2616*f3 - 1274*f4 + 251*f5);
        end

        % 计算误差
        errors(j) = abs(x(end) - x_ref);
    end

    % 计算收敛阶
    orders = log2(errors(1:end-1) ./ errors(2:end));

    % 输出误差与收敛阶
    fprintf('   N        error           order\n');
    for j = 1:length(Ns)
        if j == 1
            fprintf('%4d   %e\n', Ns(j), errors(j));
        else
            fprintf('%4d   %e     %.2f\n', Ns(j), errors(j), orders(j-1));
        end
    end

    % 可选绘图
    loglog(Ns, errors, '-o');
    grid on;
    xlabel('N');
    ylabel('误差');
    title('Adams-Bashforth 5 阶方法误差收敛');
end
