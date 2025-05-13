function query_solution(xs, ys)
    while true
        prompt = '请输入一个 x ∈ [1, %.4f] 查询对应 y 值 (输入回车退出):\n';
        xq = input(sprintf(prompt, xs(end)));

        if isempty(xq) || ~isnumeric(xq)
            break;
        end

        if xq < xs(1) || xq > xs(end)
            fprintf('超出范围 [1, %.4f].\n', xs(end));
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

        fprintf('插值结果: y(%.4f) ≈ %.10f\n', xq, yq);
    end
end
