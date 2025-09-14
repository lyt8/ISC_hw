N = 1500;
M = 3000;
x_plot = linspace(-1, 1 ,M+1);
y_plot = runge(x_plot);

tol = 1e-10;   % tolerance                   

% Compute splines for 4 cases    (Clamped / Natural) × (Uniform / Chebyshev)
% 'c' = clamped spline, 'n' = natural spline
% 'u' = uniform nodes, 'c' = Chebyshev nodes
[N_cu, all_err_cu, y_uniform_clamped] = findN(x_plot, y_plot, N, M, tol, 'clamped', 'uniform');
[N_nu, all_err_nu, y_uniform_natural] = findN(x_plot, y_plot, N, M, tol, 'natural', 'uniform');
[N_cc, all_err_cc, y_cheby_clamped] = findN(x_plot, y_plot, N, M, tol, 'clamped', 'chebyshev');
[N_nc, all_err_nc, y_cheby_natural] = findN(x_plot, y_plot, N, M, tol, 'natural', 'chebyshev');

% Plot spline curves
plot_spline(x_plot, y_plot, y_uniform_clamped, 'clamped', 'uniform', tol, N_cu)
plot_spline(x_plot, y_plot, y_uniform_natural, 'natural', 'uniform', tol, N_nu)
plot_spline(x_plot, y_plot, y_cheby_clamped, 'clamped', 'chebyshev', tol, N_cc)
plot_spline(x_plot, y_plot, y_cheby_natural, 'natural', 'chebyshev', tol, N_nc)

% Plot interpolation errors
plot_error(x_plot, all_err_cu, 'clamped', 'uniform')
plot_error(x_plot, all_err_nu, 'natural', 'chebyshev')
plot_error(x_plot, all_err_cc, 'clamped', 'uniform')
plot_error(x_plot, all_err_nc, 'natural', 'chebyshev')

% ========================================================================
function y = runge(x)
    y = 1 ./ (1 + 25 * x.^2);
end

function dy = runge_derivative(x)
    dy = -50*x ./ (1 + 25*x.^2).^2;
end

% Find the number of nodes such that the error less than tolerance
function[N, all_err, y_cubic] = findN(x_plot, y_plot, N, M, tol, bcType, nodeType)
    while true
        %uniform or chebyshev
        if strcmp(nodeType, 'uniform')
            x_nodes = linspace(-1,1,N+1);    % uniform nodes
        else
            k = 0:N;    
            x_nodes = cos((k)*pi/ N);    % chebyshev nodes
            x_nodes = sort(x_nodes);
        end
        y_nodes = runge(x_nodes);
        
        % Build and evaluate spline
        y_cubic = cubic_spline(x_plot, x_nodes, y_nodes, N, M, bcType);
        % Compute error
        all_err = y_plot - y_cubic;    
        err = max(abs(all_err)); 

        % Check tolerance
        if err < tol
            fprintf('%s spline (%s nodes) : When N = %d, the error = %e and it is less than %e', ...
                bcType, nodeType, N, err, tol);
            fprintf('\n')
            break;
        else
            N = N + 1; 
        end

        if N > M 
            fprintf('N exceeds maximum limit of %d, stopping iteration.\n', M);
            break;
        end
    end
end

% cubic spline interpolation 
function  y_cubic = cubic_spline(x_plot, x_nodes, y_nodes, N, M, bcType)
    fpo = runge_derivative(x_nodes(1));
    fpn = runge_derivative(x_nodes(end));
    h = diff(x_nodes);

    % Build system matrix & RHS
    if strcmp(bcType, 'clamped')        % clamped
        main_diag = [2*h(1), 2*(h(1:N-1) + h(2:N)), 2*h(N)];
        lower_diag = [h, 0];
        upper_diag = [0, h];

        A = spdiags([lower_diag(:), main_diag(:), upper_diag(:)], -1:1, N+1, N+1);

        rhs = [ ...
        3/h(1) * (y_nodes(2)-y_nodes(1)) - 3*fpo ; ...
        (3./h(2:N) .* (y_nodes(3:N+1)-y_nodes(2:N)) ...
        - 3./h(1:N-1) .* (y_nodes(2:N)-y_nodes(1:N-1)))' ; ...
        3*fpn - 3/h(N) * (y_nodes(N+1)-y_nodes(N)) ...
        ];
    else        % natural
        main_diag = [1, 2*(h(1:N-1) + h(2:N)), 1];
        lower_diag = [0, h];
        upper_diag = [h, 0];

        A = spdiags([lower_diag(:), main_diag(:), upper_diag(:)], -1:1, N+1, N+1);

        rhs = [ ...
        0 ; ...
        (3./h(2:N) .* (y_nodes(3:N+1)-y_nodes(2:N)) ...
        - 3./h(1:N-1) .* (y_nodes(2:N)-y_nodes(1:N-1)))' ; ...
        0 ; ...
        ];
    end

    % coefficients
    c = A \ rhs;
    b = zeros(N, 1);
    d = zeros(N, 1);
    y_cubic = zeros(1, M+1);

    for i =1:N
        b(i) = (y_nodes(i+1) - y_nodes(i)) / h(i) - h(i) * (2*c(i) + c(i+1)) / 3;
        d(i) = (c(i+1) - c(i)) / (3*h(i));
    end

    % evaluate spline
    for k = 1:length(x_plot)
        j = find(x_nodes <= x_plot(k), 1, 'last');
        if j == N+1, j = N; end
        dx = x_plot(k) - x_nodes(j);
        y_cubic(k) = y_nodes(j) + b(j)*dx + c(j)*dx^2 + d(j)*dx^3;
    end
end

% plot spline
function plot_spline(x_plot, y_plot, y_spline, bcType, nodeType, tol, N)
    figure;
    plot(x_plot, y_plot, 'LineWidth',2); hold on;
    plot(x_plot, y_spline, '--', 'LineWidth', 2);

    legend('Runge', sprintf('%s spline (%s nodes)', bcType, nodeType), 'Location', 'best');
    title(sprintf('%s cubic spline (%s nodes) with tol = %.1e, N = %d', ...
                  bcType, nodeType, tol, N));
    grid on;
end

% plot error
function plot_error(x_plot, errorData, bcType, nodeType)
    figure;
    plot(x_plot, errorData, 'LineWidth',2);
    title(sprintf('Error of %s spline (%s nodes)', bcType, nodeType))
    grid on;
end