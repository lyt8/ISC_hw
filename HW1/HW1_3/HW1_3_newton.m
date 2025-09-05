x_sin = linspace(0, 1, 2001);
y_sin = f(x_sin);

N = 10;
x_IN = linspace(0, 1, N+1);
y_IN = f(x_IN);

new_coe = divided_difference(x_IN, y_IN, N);
new_val = newton_poly(x_sin, x_IN, new_coe, N);

plot(x_sin, y_sin, 'LineWidth', 4); hold on
plot(x_sin, new_val, 'LineWidth', 2); hold on
plot(x_IN, y_IN, 'ko', 'MarkerFaceColor', 'y', 'MarkerSize', 6);
legend('f(x)', 'newton')


function y = f(x)
    y = sin(x);
end

function new_coe = divided_difference(x_IN, y_IN, N)
    F_ij = ones(N+1, N+1);
    % Compute the divided differences
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
    new_val = 0;
    prod = 1;

    for i = 1:N
        prod = prod .* (x_sin - x_IN(i));
        new_val = new_val + new_coe(i+1) * prod;
    end
    new_val = new_val + new_coe(1);
end
