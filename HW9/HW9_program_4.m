format long;
M = 10000;
h = 1/M;
plot_x = linspace(0, 1, M + 1);
plot_y = exact_u(plot_x, h);
alpha = 1 / M * (sum(f(plot_x)) - 0.5 * f(plot_x(1)) - 0.5 * f(plot_x(M+1)));
n_list = 2 .^ (3:9) ;  
errors = zeros(size(n_list));

for k = 1:length(n_list)
    n = n_list(k);
    fdm_x = linspace(0, 1, n + 1);
    fdm_u = finite_difference_method(n, @f, alpha);

    u_on_fem_x = exact_u(fdm_x);
    maxerr   = max(abs(fdm_u - (u_on_fem_x).'));
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
    y = exp(sin(x));
end

function plot_y = exact_u(x)
    plot_y = zeros(size(x));
    
    for k = 1:length(x)
        xx = x(k);
        t  = linspace(0, xx, 400);          
        f = (xx - t) .* exp(sin(t));       
        plot_y(k) = trapz(t, f);          
    end 
end

function u = finite_difference_method(n, f, alpha)
    h = 1 / n;
    x = linspace(0, 1, n+1);       
    fx = f(x);                     
    A = zeros(n, n);
    b = zeros(n, 1);

    A(1,1) = -2 / h^2;
    A(1,2) =  1 / h^2;
    b(1)   =  fx(2);        

    for i = 2:n-1
        A(i,i-1) =  1 / h^2;
        A(i,i)   = -2 / h^2;
        A(i,i+1) =  1 / h^2;
        b(i)     =  fx(i+1);   
    end

    % -----boundary-----
    A(n,n)   =  3 / (2*h);
    A(n,n-1) = -4 / (2*h);
    A(n,n-2) =  1 / (2*h);
    b(n)     =  alpha;

    u_inner = A \ b;        

    u = [0; u_inner];
end