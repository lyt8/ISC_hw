format long;
M = 4096;
plot_x = linspace(0, 1, M + 1);
dt =  1 / 10 ^ (floor(log10(2 * M^2)) + 1);
mu = dt * M ^ 2;
ref_u = mol_step(plot_x, M);

n_list = 2 .^ (3:9);
n_length = length(n_list);
error_mol = zeros(1, n_length);
error_fdm = zeros(1, n_length);

for i = 1:n_length
    n = n_list(i);
    mu = dt * n_list(i) ^ 2;
    u_x = linspace(0, 1, n_list(i) + 1);
    fdm_u = fdm_step(u_x, mu, dt, n);

    u_on_fem_x = interp1(plot_x, ref_u, u_x);
    mol_u = mol_step(u_x, n);

    error_fdm(i) = max(abs(u_on_fem_x - fdm_u));
    error_mol(i) = max(abs(u_on_fem_x - mol_u));
end

%-------------------------output-----------------------------
for i = 1: n_length
    if i == 1
        fprintf('error of fdm = %6d, order = none\n', error_fdm(1));
    else
        fprintf('error of fdm = %6d, order = %6d\n', error_fdm(i), log(error_fdm(i-1) / error_fdm(i)) / log(2));
    end
end

for i = 1: n_length
    if i == 1
        fprintf('error of mol = %6d, order = none\n', error_mol(1));
    else
        fprintf('error of mol = %6d, order = %6d\n', error_mol(i), log(error_mol(i-1) / error_mol(i)) / log(2));
    end
end

%----------figure----------

figure;
plot(plot_x, ref_u); hold on;
plot(u_x, fdm_u, '--');
legend('ref', 'fdm');
title('Exact Solution vs. Finite Difference Approximation');
grid on;

figure;
plot(plot_x, ref_u); hold on;
plot(u_x, mol_u, '--');
legend('ref', 'mol');
title('Exact Solution vs. Method of line');
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

%-------Method of line----------
function mol_u = mol_step(u_x, n)
    h = 1 / n;

    %--------v = V^-1 * u(0)-----------
    %--------use fft--------------
    mol_u_0 = ic(u_x(2:n));
    v_0 = [0, mol_u_0, 0, -flip(mol_u_0)];
    fft_mol_u_0 = (2 / n * fft(v_0)); 
    mol_u_0_hat = imag(fft_mol_u_0(2:n)) / 2;

    %-------lambda * v-----
    mol_u_hat = zeros(1, n - 1);
    for i = 1:(n-1)
        lambda_i = (-4 / h^2) * sin(i*pi/(2*n))^2;
        mol_u_hat(i) = exp(lambda_i) * mol_u_0_hat(i);
    end

    %-------V * u(t)-------------
    %--------use fft--------------
    v_1 = [0, mol_u_hat, 0, -flip(mol_u_hat)];
    fft_mol_u = fft(v_1);
    mol_u = imag(fft_mol_u(2:n)) / 2;
    mol_u = [0, mol_u, 0];
end