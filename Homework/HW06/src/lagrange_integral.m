function I = lagrange_integral(f, x_nodes, a, b)
    % 计算所有基函数的积分系数
    integrals = precompute_integrals(x_nodes, a, b);
    
    % 初始化积分值
    I = 0;
    N = length(x_nodes);
    
    % 对每个节点进行加权求和
    for i = 1:N
        % 使用预计算的基函数积分值
        I = I + f(x_nodes(i)) * integrals(i);
    end
end

function integrals = precompute_integrals(x_nodes, a, b)
    % 计算每个基函数的积分并返回
    N = length(x_nodes);
    integrals = zeros(1, N);
    
    for i = 1:N
        % 计算基函数的分母
        den = 1;
        for j = 1:N
            if j ~= i
                den = den * (x_nodes(i) - x_nodes(j));
            end
        end

        x_nodes_new = x_nodes;
        x_nodes_new(i) = [];  % 去掉第 i 个节点
        
        % 计算第 i 个基函数的积分
        integrals(i) = prod_poly_integral(x_nodes_new, a, b) / den;
    end
end

function I = prod_poly_integral(x_nodes, a, b)
    % 计算多项式 (x - x1)(x - x2)...(x - xn) 在区间 [a, b] 上的积分
    p = poly(x_nodes);  % 多项式系数，最高次在前
    N = length(p) - 1;
    I = 0;
    for k = 0:N
        coeff = p(end - k);  % x^k 的系数
        I = I + coeff * (b^(k+1) - a^(k+1)) / (k+1);
    end
end
