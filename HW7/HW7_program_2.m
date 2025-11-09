format long

n_list = [10,20,40,80,160,320,640];
num_n = length(n_list);
res_n = zeros(1,num_n);
err_n = zeros(1,num_n-1);

figure; hold on;
for k = 1:length(n_list)
    fdm_x = linspace(0, 1, n_list(k)+1);
    fdm_u_0 = ones(1, n_list(k)+1);
    res = zeros(1, 6);

    [fdm_res, fdm_u] = finite_difference_method(fdm_u_0, n_list(k));
    res(1) = fdm_res;
    for i = 1:5
        [fdm_res, fdm_u] = finite_difference_method(fdm_u, n_list(k));
        res(i) = fdm_res;
    end

    res_n(k) = fdm_res;
    if k >= 2
        fdm_x_coarse = linspace(0, 1, n_list(k-1)+1);
        u_fine_interp = interp1(fdm_x, fdm_u, fdm_x_coarse, 'linear');
        err = max(abs(u_fine_interp - fdm_u_old));
        err_n(k-1) = err ;
    end
    fdm_u_old = fdm_u;
    plot(fdm_x, fdm_u); hold on;

    fprintf('n = %4d, res = %.6e\n', n_list(k), fdm_res);
end
grid on;
xlabel('x'); 
ylabel('u');

% Calculate and display the error norms
    for k = 2: num_n -1
        fprintf('Error between n = %4d and n = %4d: %.6e\n', n_list(k-1), n_list(k), err);
    end

%----------function-------------
function [norm_res,fdm_u] = finite_difference_method(u_0, n) 
    h = 1 / n; 
    fdm_u = u_0;
    res = ((-fdm_u(1:n-1) + 2*fdm_u(2:n) - fdm_u(3:n+1)) / h^2 + sin(fdm_u(2:n)));
    norm_res = max(abs(res)) ;

    %------A = LU------
    L = zeros(1, n-1);
    U = zeros(1, n-2);
    L(1) = 2 + h^2 * cos(u_0(1));

    for i = 1:n-2
        U(i) = -1 / L(i);
        L(i+1) = 2 + h^2 * cos(u_0(i+1)) + U(i); 
    end

    %------Lz = f------
    z = zeros(n-1, 1); 
    z(1) = -res(1) * h^2 / L(1);
    for i = 2:n-1
        z(i) = (-res(i) * h^2 + z(i-1)) / L(i);
    end

    %------Ux = z------
    x = zeros(1,n-1); 
    x(n-1) = z(n-1);
    for i = 2:n-1
        x(n - i) = z(n-i) - U(n - i) * x(n + 1 - i);
    end

    detla = [0, x, 0];
    fdm_u = fdm_u + detla;


end