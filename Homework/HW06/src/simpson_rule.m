function I_simpson = simpson_rule(f, a, b, N)
    h = (b - a) / N;
    x = a:h:b;  % 节点
    fx = f(x);
    I_simpson = (h / 3) * (fx(1) + fx(end) + 4 * sum(fx(2:2:end)) + 2 * sum(fx(3:2:end-1)));
end
