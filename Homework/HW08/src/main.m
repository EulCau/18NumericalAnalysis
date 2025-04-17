disp('选择函数进行 Romberg 积分: ');
disp('1. sin(x)/x from 0 to 1');
disp('2. (cos(x)-exp(x))/sin(x) from -1 to 1');
disp('3. 1/(x*exp(x)) from 1 to ∞');

choice = input('输入 1, 2, 或 3: ');

switch choice
    case 1
        f = @(x) f1_with_limit(x);
        a = 0; b = 1;
    case 2
        f = @(x) f2_with_limit(x);
        a = -1; b = 1;
    case 3
        f = @(x) f3_with_limit(x);
        a = 0; b = 1;
end

max_iter = input('输入计算的行数: ');

R = romberg(f, a, b, max_iter);

% 显示结果
fprintf('\nRomberg 表格（前 %d 行）：\n', max_iter);
for i = 1:max_iter
    for j = 1:i
        fprintf('%15.10f ', R(i,j));
    end
    fprintf('\n');
end

% --- 特殊处理 sin(x)/x 函数 ---
function y = f1_with_limit(x)
    y = ones(size(x));  % 限值为 1
    near_zero = abs(x) < 1e-8;
    y(~near_zero) = sin(x(~near_zero)) ./ x(~near_zero);
end

% --- 特殊处理 (cos(x) - e^x)/sin(x) 函数 ---
function y = f2_with_limit(x)
    y = -ones(size(x));  % 限值为 -1
    near_zero = abs(x) < 1e-8;
    y(~near_zero) = (cos(x(~near_zero)) - exp(x(~near_zero))) ./ sin(x(~near_zero));
end

% --- 特殊处理 e^(-1/t)/t 函数 ---
function y = f3_with_limit(x)
    y = zeros(size(x));  % 限值为 0
    near_zero = abs(x) < 1e-8;
    y(~near_zero) = exp(-1./x(~near_zero))./x(~near_zero);
end
