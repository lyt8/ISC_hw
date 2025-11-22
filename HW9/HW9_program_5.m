format long;
eps = 0.01;
M = 10000;
plot_x = linspace(0, 1, M + 1);
plot_y = exact_u(plot_x, eps);

n_list = 2 .^ (6:21) ;  
errors = zeros(size(n_list));

for k = 1:length(n_list)
    n = n_list(k);
    fdm_x = linspace(0, 1, n + 1);
    fdm_u = finite_difference_method(n, eps);

    u_on_fem_x = exact_u(fdm_x, eps);
    maxerr   = max(abs(fdm_u - u_on_fem_x));
    errors(k) = maxerr;
    fprintf('n = %6d, h = %.3e, max error = %.6e\n', n, 1/n, errors(k));
end

%--------------figure--------------

figure;
plot(plot_x, plot_y, 'LineWidth',3); hold on;
plot(fdm_x, fdm_u, '--', 'LineWidth',3);
legend('Reference Solution', 'Finite Difference Method');
xlabel('x');
ylabel('u(x)');
title('Exact Solution vs. Finite Difference Approximation');
grid on;

figure;
loglog(1./n_list, errors);
xlabel('h = 1/n  (grid size)');
ylabel('Maximum Error  ||u_h - u_{exact}||_\infty');
title('Error vs. Mesh Size (Finite Difference Method)');
grid on;

%-------------function--------------
function y = f(x)
    y = zeros(1, lenght(x));
end

function plot_y = exact_u(x, eps)
    plot_y = (-exp(1) / (1 - exp(-1 / eps + 1))) .* exp(-1 / eps .* x) + ...
        (exp(1) / (1 - exp(-1 / eps + 1)) .* exp(-1 .* x)); 
end

function fdm_u = finite_difference_method(n, eps) 
    h = 1 / n;
    A_tri = h^2 - 2 * eps;
    A_lower = eps - h / 2 - h * eps / 2;
    A_upper = eps + h / 2 + h * eps / 2;
    rhs = [0, -eps - h / 2 - h * eps / 2];

    %------A = LU------
    L = zeros(1, n-2);
    U = zeros(1, n-1);
    U(1) = A_tri;

    for i = 1:n-2
        L(i) = A_lower / U(i);
        U(i + 1) = A_tri - L(i) * A_upper;
    end

    %------Lz = f------
    z = zeros(n - 1, 1); 
    z(1) = rhs(1);

    for i = 3:n-1
        z(i - 1) = rhs(1) - z(i - 2) * L(i - 2);
    end
    z(n - 1) = rhs(2) - z(n - 2) * L(n - 2);
    
    %------Ux = z------
    x = zeros(1, n - 1); 
    x(n - 1) = z(n - 1) / U(n - 1);
    
    for i = 2:n-2
        x(n - i) = (z(n - i) - A_upper * x(n - i + 1)) / U(n - i);
    end
    x(1) = (z(1) - A_upper * x(2)) / U(1);
    
    fdm_u = [0, x, 1];
end