function rkf_adaptive_solver()
    % 初始条件与控制参数
    x0 = 1;
    y0 = 1;
    h = 0.01;
    tol = 1e-5;
    x = x0;
    y = y0;

    % 存储解轨迹
    xs = x;
    ys = y;

    % 定义微分方程
    f = @(x, y) exp(y*x) - cos(y - x);

    k = 0;
    while k <= 1000000
        k = k+1;
        try
            % 计算 RKF45 的 6 阶和 5 阶近似
            [y4, y5] = rkf45_step(f, x, y, h);

            % 局部误差估计
            err = abs(y5 - y4);

            % 步长控制（策略2）
            if err < tol
                % 接受该步
                x = x + h;
                y = y5;  % 用高阶值更新

                % 存储
                xs(end+1) = x;
                ys(end+1) = y;
            end

            % 计算新步长
            delta = (tol / max(err, 1e-14))^(1/5);  % p=4
            h = 0.9 * h * min(5, max(0.2, delta));

            if ~isfinite(y)
                disp('检测到数值溢出，终止计算。');
                break;
            end
        catch
            disp('数值溢出或函数错误，终止。');
            break;
        end
    end

    % 输出解范围
    fprintf('解计算完成。\n');
    fprintf('解的范围为 x ∈ [%.4f, %.4f]\n', xs(1), xs(end));

    % 用户插值查询
    while true
        prompt = '请输入一个 x ∈ [1, %.4f] 查询对应 y 值（输入非数字退出）:\n';
        xq = input(sprintf(prompt, xs(end)));

        if isempty(xq) || ~isnumeric(xq)
            break;
        end

        if xq < xs(1) || xq > xs(end)
            fprintf('超出范围 [1, %.4f]。\n', xs(end));
            continue;
        end

        % 找邻近的两点进行线性插值
        idx = find(xs <= xq, 1, 'last');
        if idx == length(xs)
            idx = idx - 1;
        end
        x1 = xs(idx); x2 = xs(idx+1);
        y1 = ys(idx); y2 = ys(idx+1);
        yq = y1 + (xq - x1) * (y2 - y1) / (x2 - x1);

        fprintf('插值结果：y(%.4f) ≈ %.10f\n', xq, yq);
    end
end

function [y4, y5] = rkf45_step(f, x, y, h)
    % Runge-Kutta-Fehlberg RKF45 经典系数
    k1 = f(x, y);
    k2 = f(x + h/4, y + h*k1/4);
    k3 = f(x + 3*h/8, y + h*(3*k1 + 9*k2)/32);
    k4 = f(x + 12*h/13, y + h*(1932*k1 - 7200*k2 + 7296*k3)/2197);
    k5 = f(x + h, y + h*(439*k1/216 - 8*k2 + 3680*k3/513 - 845*k4/4104));
    k6 = f(x + h/2, y + h*(-8*k1/27 + 2*k2 - 3544*k3/2565 + 1859*k4/4104 - 11*k5/40));

    % 四阶近似（低精度）
    y4 = y + h*(25*k1/216 + 1408*k3/2565 + 2197*k4/4104 - k5/5);

    % 五阶近似（高精度）
    y5 = y + h*(16*k1/135 + 6656*k3/12825 + 28561*k4/56430 - 9*k5/50 + 2*k6/55);
end
