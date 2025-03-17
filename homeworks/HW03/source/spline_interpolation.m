function spline_interpolation()
    N_vals = [5, 10, 20, 40];
    errors_linear = zeros(size(N_vals));
    errors_cubic = zeros(size(N_vals));
    
    for j = 1:length(N_vals)
        N = N_vals(j);
        x = linspace(0, 1, N+1); 
        y = exp(x);
        x_mid = (x(1:end-1) + x(2:end)) / 2;
        
        % Linear spline interpolation
        S1 = linear_spline(x, y, x_mid);
        errors_linear(j) = max(abs(exp(x_mid) - S1));
        
        % Cubic spline interpolation with given boundary conditions
        S3 = cubic_spline(x, y, x_mid);
        errors_cubic(j) = max(abs(exp(x_mid) - S3));
    end
    
    % Compute order of convergence
    ord_linear = compute_order(errors_linear, N_vals);
    ord_cubic = compute_order(errors_cubic, N_vals);
    
    % Display results
    disp('N    Error(Linear)   Order(Linear)   Error(Cubic)   Order(Cubic)');
    for i = 1:length(N_vals)
        fprintf('%d    %.5e    %.2f    %.5e    %.2f\n', N_vals(i), errors_linear(i), ord_linear(i), errors_cubic(i), ord_cubic(i));
    end
end

function S = linear_spline(x, y, x_mid)
    S = zeros(size(x_mid));
    for i = 1:length(x)-1
        idx = (x_mid >= x(i)) & (x_mid <= x(i+1));
        S(idx) = y(i) + (y(i+1) - y(i)) / (x(i+1) - x(i)) * (x_mid(idx) - x(i));
    end
end

function S = cubic_spline(x, y, x_mid)
    N = length(x) - 1;
    h = diff(x);

    A = 2 * diag(ones(N-1,1),0) + diag(h(1:end-1)./ (h(1:end-1) + h(2:end)), -1) + diag(h(2:end)./ (h(1:end-1) + h(2:end)), 1);
    rhs = 6 * diff((diff(y) ./ h) ./ h);
    
    % Apply boundary conditions
    d0 = 0;  % S'(x_0) = 0
    dN = exp(1);  % S'(x_N) = e
    
    M = A \ rhs;  % Solve for second derivatives
    M = [d0; M; dN];
    
    % Compute cubic spline values
    S = zeros(size(x_mid));
    for i = 1:N
        idx = (x_mid >= x(i)) & (x_mid <= x(i+1));
        xi = x(i);
        xi1 = x(i+1);
        hi = xi1 - xi;
        
        S(idx) = ((M(i) * (xi1 - x_mid(idx)).^3 + M(i+1) * (x_mid(idx) - xi).^3) / (6*hi)) ...
               + ((y(i)/hi - M(i)*hi/6) * (xi1 - x_mid(idx))) ...
               + ((y(i+1)/hi - M(i+1)*hi/6) * (x_mid(idx) - xi));
    end
end

function ord = compute_order(errors, N_vals)
    ord = zeros(size(errors));
    for i = 2:length(N_vals)
        ord(i) = log(errors(i-1) / errors(i)) / log(N_vals(i) / N_vals(i-1));
    end
end
