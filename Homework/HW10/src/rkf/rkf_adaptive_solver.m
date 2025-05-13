function rkf_adaptive_solver(f, x0, y0, ifplot)
    if isempty(ifplot)
        ifplot = true;
    end

    % 初始控制参数
    h = 0.01;
    tol = 1e-5;
    N = 65536;
    k = 2;

    % 存储解轨迹
    xs = zeros(1, N);
    ys = zeros(1, N);
    xs(1) = x0;
    ys(1) = y0;

    while k <= N
        try
            x = xs(k - 1);
            y = ys(k - 1);
            % 计算 RKF45 的 6 阶和 5 阶近似
            [y4, y5] = rkf45_step(f, x, y, h);

            % 局部误差估计
            err = abs(y5 - y4);

            % 步长控制 (策略 2)
            if err < tol
                % 接受该步
                x_old = x;
                x = x + h;
                if x == x_old
                    disp('检测到x无法累加, 终止计算.');
                    break;
                end
                y = y5;  % 用高阶值更新

                % 存储
                xs(k) = x;
                ys(k) = y;
                k  = k + 1;
            end

            % 计算新步长
            delta = (tol / max(err, 1e-14))^(1/5);  % p=4
            h = 0.9 * h * min(5, max(0.2, delta));

            if ~isfinite(y)
                disp('检测到数值溢出, 终止计算.');
                break;
            end
        catch
            disp('数值溢出或函数错误, 终止.');
            break;
        end
    end
    if any(ys ~= 0)
        idx = find(ys ~= 0, 1, 'last');
        xs = xs(1:idx);
        ys = ys(1:idx);
    else
        xs = x0;
        ys = y0;
    end

    % 输出解范围
    fprintf('解计算完成.\n');
    fprintf('解的范围为 x ∈ [%.4f, %.32f]\n', xs(1), xs(end));
    if ifplot
        plot(xs, ys);
        grid on;
    end

    % 用户插值查询
    query_solution(xs, ys)
end
