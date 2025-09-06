x_sin = linspace(0, 1, 2001);
y_sin = f(x_sin);

N = 10;
x_IN = linspace(0, 1, N+1);
y_IN = f(x_IN);

stand_val = stand_poly(x_sin, x_IN, y_IN, N);

new_coe = divided_difference(x_IN, y_IN, N);
new_val = newton_poly(x_sin, x_IN, new_coe, N);

lan_val = lagrange_basis_func(x_sin, x_IN, y_IN, N);

bar_val = barycentric_lagrange_poly(x_sin, x_IN, y_IN, N);

plot(x_sin, y_sin, 'LineWidth', 4); hold on
plot(x_sin, stand_val, 'LineWidth', 2); hold on
plot(x_sin, new_val, 'LineWidth', 2); hold on
plot(x_sin, lan_val, 'LineWidth', 2); hold on
plot(x_sin, bar_val, 'LineWidth', 2); hold on
plot(x_IN, y_IN, 'ko', 'MarkerFaceColor', 'y', 'MarkerSize', 6);

xlabel('x'); ylabel('y'); 
legend('f(x)', 'standard', 'newton', 'lagrange', 'barycentric')
grid on

function y = f(x)
    y = sin(x);
end

% standard basis
function stand_val = stand_poly(x_sin, x_IN, y_IN, N)
    stand_val = zeros(size(x_sin));
    ven_matrix = ones(N+1, N+1);
    for i = 1:N+1
        for j = 1:N+1
            ven_matrix(i, j) = x_IN(i) ^ (j-1);
        end
    end

    coe = ven_matrix \ y_IN(:);
    for j = 1:N+1
        stand_val = stand_val + coe(j) * x_sin .^ (j-1);
    end
end

%Newton interpolation
function new_coe = divided_difference(x_IN, y_IN, N)
    F_ij = ones(N+1, N+1);
    for j = 1:N+1
        F_ij(j, 1) = y_IN(j);
    end
    for j = 2:N+1
        for i = 1:N+2-j
            F_ij(i, j) = (F_ij(i+1, j-1) - F_ij(i, j-1)) / (x_IN(i+j-1) - x_IN(i));
        end
    end
    new_coe = F_ij(1,:);
end

function new_val = newton_poly(x_sin, x_IN, new_coe, N)
    new_val = new_coe(N+1);
    for k = N:-1:1
        new_val = new_coe(k) + (x_sin - x_IN(k)) .* new_val;
    end
end

%lagrange interpolation
function lan_val = lagrange_basis_func(x_sin, x_IN, y_IN, N)
    Li = ones(N, 2001); 
    for i = 1:N
        for j = [1:i-1, i+1:N]
            Li(i,:) = Li(i,:) .*  (x_sin - x_IN(j)) / (x_IN(i) - x_IN(j));
        end
        Li(i,:) = Li(i,:) * y_IN(i);
    end
    lan_val = sum(Li, 1);
end

%barycentric lagrange interpolation
function bar_val = barycentric_lagrange_poly(x_sin, x_IN, y_IN, N)
    bar_val = ones(1, 2001);
    wi = ones(1, N+1);
    for j = 1:N+1
        wi(j) = wi(j) * prod(x_IN(j) - x_IN([1:j-1, j+1:N+1])); 
    end
    wi = 1 ./ wi;

    for i = 1:2001
        [dmin, j] = min(abs(x_sin(i) - x_IN));
        if dmin < 1e-6
            bar_val(i) = y_IN(j);      
        else
            denom = sum( wi(1:N+1) ./ (x_sin(i) - x_IN(1:N+1)) );
            numer = sum( (wi(1:N+1) .* y_IN(1:N+1)) ./ (x_sin(i) - x_IN(1:N+1)) );
            bar_val(i) = numer / denom;
        end
    end
end