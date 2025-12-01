format long;
M = 4096;
plot_x = linspace(0, 1, M + 1);
dt =  1 / M^2 / 2;
mu = dt * M ^ 2;
ref_u = fdm_step(plot_x, mu, dt, M);

t_list = 1 ./ 2 .^ (10:15);
t_length = length(t_list);
error_fdm = zeros(1, t_length);
n = 8;

for i = 1:t_length
    dt = t_list(i);
    mu = dt * n ^ 2;
    u_x = linspace(0, 1, n + 1);
    fdm_u = fdm_step(u_x, mu, dt, n);

    u_on_fem_x = interp1(plot_x, ref_u, u_x);

    error_fdm(i) = max(abs(u_on_fem_x - fdm_u));
end

%-------------------------output-----------------------------
for i = 1: t_length
    if i == 1
        fprintf('error of fdm = %6d, order = none\n', error_fdm(1));
    else
        % 2.073032e-06 is error of fdm at n = 8 by HW10_program_2
        e_b = abs((error_fdm(i-1) - (2.073032e-06)));
        e_a = abs((error_fdm(i) - (2.073032e-06)));
        fprintf('error of fdm_t = %6d, order = %6d\n', e_a, log(e_b / e_a) / log(2));
    end
end

%----------figure----------

figure;
plot(plot_x, ref_u); hold on;
plot(u_x, fdm_u, '--');
legend('ref', 'fdm');
title('Exact Solution vs. Finite Difference Approximation');
grid on;

%---------function-----------
function u_ic = ic(x)
    u_ic = sin(2 * pi .* x) .* exp(x);
end

%----------finite diffeerence method------------
function fdm_u =  fdm_step(u_x, mu, dt, n)
    t_times = 1 / dt;

    %--------v = V^-1 * u(0)-----------
    %--------use fft--------------
    fdm_u_0 = ic(u_x(2:n));
    v_0 = [0, fdm_u_0, 0, -flip(fdm_u_0)];
    fft_fdm_u_0 = (2 / n * fft(v_0)); 
    fdm_u_0_hat = imag(fft_fdm_u_0(2:n)) / 2;

    %-------lambda * v-----
    fdm_u_hat = zeros(1, n - 1);
    for i = 1:(n-1)
        lambda = 1 - 4 * mu * sin(i*pi/(2*n))^2;
        fdm_u_hat(i) = lambda ^ t_times * fdm_u_0_hat(i);
    end

    %-------V * u(t)-------------
    %--------use fft--------------
    v_1 = [0, fdm_u_hat, 0, -flip(fdm_u_hat)];
    fft_u = fft(v_1);
    fdm_u = imag(fft_u(2:n)) / 2;
    fdm_u = [0, fdm_u, 0];
    
end

% error of fdm = 1.269704e-06, order = none
% error of fdm_t = 4.050533e-07, order = 9.878776e-01
% error of fdm_t = 2.033664e-07, order = 9.940306e-01
% error of fdm_t = 1.018819e-07, order = 9.971838e-01
% error of fdm_t = 5.097890e-08, order = 9.989252e-01
% error of fdm_t = 2.548718e-08, order = 1.000129e+00