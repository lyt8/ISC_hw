format long;
M = 100;
n = 32;
mu = [0.5, 0.509];

t_times = fix(1 / (mu(1) * (1 / M) ^ 2));
plot_x = linspace(0, 1, M + 1);
plot_t = linspace(0, 1, t_times);      

fdm_x = linspace(0, 1, n + 1);
[X, T] = meshgrid(plot_x, plot_t);
U = f(X, T);                  

[error_1, fdm_u_1, u_on_fdm_x] = fdm_step(mu(1), n);
[error_2, fdm_u_2, u_on_fdm_x] = fdm_step(mu(2), n);
fprintf('error of mu = 0.5 : %6d, error of mu = 0.509 : %6d\n', error_1, error_2);

%----------figure-------------
draw_re_fdm(fdm_x, fdm_u_1, u_on_fdm_x);
draw_re_fdm(fdm_x, fdm_u_2, u_on_fdm_x);

figure;
surf(X,T,U);
shading interp;
xlabel('x');
ylabel('t');
zlabel('u(x,t)');
title('Exact Solution Surface Plot');
colorbar;
view(135,30);

%-------------function--------------
function plot_y = f(x,t)
    plot_y = exp(-pi^2 / 4 * t ) .* sin(0.5 * pi .* x) + ...
        exp(-4 * pi^2 * t) * 0.5 .* sin(2 * pi .* x);
end

function fdm_u_ic = ic(x)
    fdm_u_ic = sin(0.5 *  pi .* x) + 0.5 * sin(2 * pi .* x);
end

function fdm_u_left = bc_left(t)
    fdm_u_left = 0;
end

function fdm_u_right = bc_right(t)
    fdm_u_right = exp(-pi^2 * t / 4);
end

function fdm_u_new = forward_euler(n, u_old, bc_left, bc_right, mu, t_times)
    fdm_u_new =zeros(1, n+1);
    t = t_times * mu * (1 / n) ^ 2; 
    A_tri = 1 - 2 * mu;
    A_lower = mu;
    A_upper = mu;

    fdm_u_new(1) = bc_left(t);
    fdm_u_new(n + 1) = bc_right(t);
    for i = 2 : n
        fdm_u_new(i) = A_tri * u_old(i) + A_lower * u_old(i - 1) + A_upper * u_old(i + 1);
    end
end

function[error, fdm_u, u_on_fdm_x] = fdm_step(mu, n)
    dt = mu * (1 / n) ^ 2;
    fdm_x = linspace(0, 1, n + 1);
    fdm_u_old = ic(fdm_x);

    t_times = fix(1 / dt);
    error_t = zeros(1, t_times);
    for j = 1:t_times
        fdm_u_new = forward_euler(n, fdm_u_old, @bc_left, @bc_right, mu, j);
        u_on_fdm_x = f(fdm_x, j * dt);
        error_t(j) = max(abs(fdm_u_new - u_on_fdm_x));

        fdm_u_old = fdm_u_new;
    end

    fdm_u = fdm_u_new;
    error = error_t(t_times);
    
end

function draw_re_fdm(fdm_x, fdm_u, u_on_fdm_x)
    figure;
    plot(fdm_x, u_on_fdm_x, 'LineWidth', 4);hold on;
    plot(fdm_x, fdm_u, '--', 'LineWidth', 4); 
    legend('Reference Solution', 'Finite Difference Method');
    xlabel('x');
    ylabel('u(x)');
    title('Exact Solution vs. Finite Difference Approximation');
    grid on;
end