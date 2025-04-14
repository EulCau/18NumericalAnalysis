function I_trap = trapezoidal_rule(f, a, b, N)
    h = (b - a) / N;
    x = a:h:b;  % 节点
    fx = f(x);
    I_trap = h * (sum(fx) - 0.5 * (fx(1) + fx(end)));  % 梯形法
end
