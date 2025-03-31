function fourier()
    N_values = [4, 8, 16, 32, 64];
    x = linspace(0, pi/2, 1000); % 计算误差时使用的高精度采样点
    
    u = @(x) 1 ./ (5 - 4 * cos(x));
    v = @(x) sin(x / 2);
    
    for N = N_values
        % 计算傅里叶系数
        u_hat = compute_fourier_coefficients(u, N);
        v_hat = compute_fourier_coefficients(v, N);
        
        % 计算离散傅里叶变换系数
        u_tilde = compute_dft_coefficients(u, N);
        v_tilde = compute_dft_coefficients(v, N);
        
        % 重建函数
        P_N_u = reconstruct_function(u_hat, x, N);
        P_N_v = reconstruct_function(v_hat, x, N);
        I_N_u = reconstruct_function(u_tilde, x, N);
        I_N_v = reconstruct_function(v_tilde, x, N);
        
        % 计算误差
        err_P_u = abs(u(x) - P_N_u);
        err_P_v = abs(v(x) - P_N_v);
        err_I_u = abs(u(x) - I_N_u);
        err_I_v = abs(v(x) - I_N_v);
        
        % 画误差图
        plot_error(x, err_P_u, err_P_v, err_I_u, err_I_v, N);
    end
end

function coeffs = compute_fourier_coefficients(f, N)
    coeffs = zeros(1, N);
    for n = -N/2:N/2
        integral_fun = @(x) f(x) .* exp(-1i * n * x);
        coeffs(n + N/2 + 1) = (1 / (2 * pi)) * integral(integral_fun, 0, 2 * pi);
    end
end

function coeffs = compute_dft_coefficients(f, N)
    x_j = linspace(0, 2 * pi, N + 1);
    x_j(end) = []; % 去掉最后一个点，防止重复
    coeffs = zeros(1, N);
    for n = -N/2:N/2
        C = 1 + (abs(n) == N/2);
        coeffs(n + N/2 + 1) = (1 / (C * N)) * sum(f(x_j) .* exp(-1i * n * x_j));
    end
end

function f_reconstructed = reconstruct_function(coeffs, x, N)
    f_reconstructed = zeros(size(x));
    for n = -N/2:N/2
        f_reconstructed = f_reconstructed + coeffs(n + N/2 + 1) * exp(1i * n * x);
    end
    f_reconstructed = real(f_reconstructed);
end

function plot_error(x, err_P_u, err_P_v, err_I_u, err_I_v, N)
    figure;
    subplot(2,1,1);
    plot(x, err_P_u, 'r', x, err_I_u, 'b');
    title(['误差 (u) N = ', num2str(N)]);
    legend('P_N u', 'I_N u');
    
    subplot(2,1,2);
    plot(x, err_P_v, 'r', x, err_I_v, 'b');
    title(['误差 (v) N = ', num2str(N)]);
    legend('P_N v', 'I_N v');
end
