format long;
M = 10000;
plot_x = linspace(0, 1, M+1);
plot_y = exact_u(plot_x);

n_list = 2 .^ (4:18);  
errors = zeros(size(n_list));  

for k = 1:length(n_list)
    n = n_list(k);
    fdm_x = linspace(0, 1, n+1);
    fdm_u = finite_difference_method(n);

    u_on_fem_x = exact_u(fdm_x);
    errors(k)   = max(abs(fdm_u - u_on_fem_x));
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
function plot_y = exact_u(x)
    plot_y = ((1 + exp(1)) / 2 / exp(1)) .* x .* exp(x)-...
        exp(x) + 1;
end

function fdm_u = finite_difference_method(n) 
    h = 1 / n; 
    A_tri = [-2 + h^2, (-2 - 4 * h + 3 * h^2) / 3];
    A_lower = [1 + h, (2 + 4 * h) / 3];
    A_upper = 1 - h;
    rhs = [h^2, (5 * h^2 - 2 * h) / 3];

    %------A = LU------
    L = zeros(1, n-2);
    U = zeros(1, n-1);
    U(1) = A_tri(1);

    for i = 1:n-3
        L(i) = A_lower(1) / U(i);
        U(i+1) = A_tri(1) - L(i) * A_upper;
    end
    L(n-2) = A_lower(2) / U(n-2);
    U(n-1) = A_tri(2) - L(n-2) * A_upper;

    %------Lz = f------
    z = zeros(n-1, 1); 
    z(1) = rhs(1);
    for i = 2:n-2
        z(i) = rhs(1) - z(i-1) * L(i-1);
    end
    z(n-1) = rhs(2) - z(n-2) * L(n-2);

    %------Ux = z------
    x = zeros(1, n-1); 
    x(n-1) = z(n-1) / U(n-1);
    for i = 2:n-1
        x(n - i) = (z(n - i) - A_upper * x(n - i + 1)) / U(n - i);
    end
    x_N = (2 * h - x(n - 2) + 4 * x(n - 1)) / 3;

    fdm_u = [0, x, x_N];
end