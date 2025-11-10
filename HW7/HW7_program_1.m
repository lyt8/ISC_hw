format long
N = 500000;
plot_x = linspace(0, 1, N+1);
y = u(plot_x, @f);

n_list = [750,1500,3000,6000,12000,24000,48000,96000,192000];  
errors = zeros(size(n_list));  

fdm_x = linspace(0, 1, n_list(1)+1);
fdm_u = finite_difference_method(fdm_x, n_list(1), @f);
rich = zeros(1, length(n_list));

for k = 1:length(n_list)
    n = n_list(k);
    fdm_x_new = linspace(0, 1, n+1);
    fdm_u_new = finite_difference_method(fdm_x_new, n, @f);

    u_on_fem_x = interp1(plot_x, y, fdm_x_new);
    maxerr   = norm(fdm_u_new - u_on_fem_x);
    errors(k) = maxerr;

    rich(k) = richardson(fdm_x_new, fdm_u_new, fdm_u, fdm_x);
    fdm_x = fdm_x_new;
    fdm_u = fdm_u_new;
    fprintf('n = %6d, h = %.3e, max error = %.6e, richard = %.4e\n', n, 1/n, maxerr,  rich(k));
end

%-----------figure-------------
figure;
plot(plot_x, y, 'LineWidth',3); hold on;
plot(fdm_x, fdm_u, '--', 'LineWidth',3);
legend('Reference Solution', 'Finite Difference Method(n=192000)');
xlabel('x');
ylabel('u(x)');
grid on;

figure;
plot(1./n_list(2:end), errors(2:end), 'o-', 'LineWidth', 3);
xlabel('Grid Spacing (h)');
ylabel('Max Error');
title('Error Analysis of Finite Difference Method');
grid on;

%-----------function-------------
function y = f(x)
    y = exp(sin(x));
end

function y = u(x, f)
    f_ref = f(x);

    I1 = cumtrapz(x, x .* f_ref);      
    I2 = cumtrapz(x, (1 - x) .* f_ref); 
    I2_1 = I2(end);

    y = (1 - x) .* I1 + x .* (I2_1 - I2);
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
        z(i) = (f(fdm_x(i+1)) * h^2 + z(i-1)) / L(i);
    end

    %------Ux = z------
    x = zeros(1,n-1); 
    x(n-1) = z(n-1);
    for i = 2:n-1
        x(n - i) = z(n-i) - U(n - i) * x(n + 1 - i);
    end

    fdm_u = [0, x, 0];
end

function rich = richardson(fdm_x_new, fdm_u_new, fdm_u, fdm_x)
    u_fine_interp = interp1(fdm_x_new, fdm_u_new, fdm_x);
    rich = norm(u_fine_interp - fdm_u) / 3;
end
