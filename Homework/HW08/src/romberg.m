% --- Romberg 算法实现 ---
function R = romberg(f, a, b, n)
    R = zeros(n, n);
    for k = 1:n
        if k == 1
            R(1,1) = (b - a)/2 * (f(a) + f(b));
        else
            h = (b - a) / 2^(k-1);
            x = a + h : 2*h : b - h;
            R(k,1) = (1/2)*R(k-1,1) + h * sum(f(x));
            for j = 2:k
                R(k,j) = R(k,j-1) + (R(k,j-1) - R(k-1,j-1)) / (4^(j-1) - 1);
            end
        end
    end
end
