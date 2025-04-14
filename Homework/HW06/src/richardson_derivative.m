function D = richardson_derivative(f, x, h, M)
    % 输入:
    % f: 被求导的函数
    % x: 计算导数的点
    % h: 步长 (初始值为1)
    % M: 外推法的次数
    % 输出:
    % D: 外推三角阵列
    
    % 初始化三角阵列 D
    D = zeros(M+1, M+1);
    
    % 计算第 0 列的差分值
    for i = 0:M
        D(i+1, 1) = forward_difference(f, x, h);
        h = h / 2;
    end
    
    % Richardson 外推计算其它列的值
    for j = 1:M
        for i = j:M
            D(i+1, j+1) = richardson_extrapolation(D(i, j), D(i+1, j), j);
        end
    end
end

% 前向差分法计算导数
function D = forward_difference(f, x, h)
    % 前向差分法计算第 i 次的导数
    D = (f(x + h) - f(x)) / h;
end

% Richardson 外推公式
function D = richardson_extrapolation(D1, D2, j)
    % Richardson 外推法
    D = (4^j * D2 - D1) / (4^j - 1);
end
