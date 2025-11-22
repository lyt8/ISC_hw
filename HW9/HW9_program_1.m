format long;
M = 10000;
plot_x = linspace(0, 1, M+1);
plot_y = exact_u(plot_x);

n_list = 4 .^ (2:12);  
errors = zeros(size(n_list));

for k = 1:length(n_list)
    n = n_list(k);
    fdm_x = linspace(0, 1, n+1);
    fdm_u = finite_difference_method(fdm_x, n, @f);

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
function y = f(x)
    y = x >= 0.4 & x <= 0.6;
end

function plot_y = exact_u(x)
    plot_y = (x < 0.4)* -0.1 .* x + ...
        (x >= 0.4 & x <= 0.6) .*(0.5 .* x.^2 -0.5 .* x + 0.08)+...
        (x > 0.6) * -0.1 .* (1 - x);
end

function fdm_u = finite_difference_method(fdm_x, n, f) 
    h = 1 / n; 
    
    %------A = LU------
    L = zeros(1, n-1);
    U = zeros(1, n-2);
    L(1) = 2;

    for i = 1:n-2
        U(i) = -1 / L(i);
        L(i+1) = 2 + U(i); 
    end

    %------Lz = f------
    z = zeros(n-1, 1); 
    z(1) = f(fdm_x(2)) * h^2 / L(1);
    for i = 2:n-1
        z(i) = (-f(fdm_x(i+1)) * h^2 + z(i-1)) / L(i);
    end

    %------Ux = z------
    x = zeros(1,n-1); 
    x(n-1) = z(n-1);
    for i = 2:n-1
        x(n - i) = z(n-i) - U(n - i) * x(n + 1 - i);
    end

    fdm_u = [0, x, 0];
end