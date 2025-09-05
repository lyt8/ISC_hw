x_sin = linspace(0, 1, 2001);
y_sin = f(x_sin);

N = 10;
x_IN = linspace(0, 1, N+1);
y_IN = f(x_IN);

stand_val = stand_poly(x_sin, x_IN, y_IN, N);

plot(x_sin, y_sin, 'LineWidth', 4); hold on
plot(x_sin, stand_val, 'LineWidth', 2); hold on
plot(x_IN, y_IN, 'ko', 'MarkerFaceColor', 'y', 'MarkerSize', 6);
legend('f(x)', 'standard', 'newton')

function y = f(x)
    y = sin(x);
end

function stand_val = stand_poly(x_sin, x_IN, y_IN, N)
    stand_val = zeros(size(x_sin));
    ven_matrix = ones(N+1, N+1);
    for i = 1:N+1
        for j = 1:N+1
            ven_matrix(i, j) = x_IN(i) ^ (j-1);
        end
    end

    coe = ven_matrix \ rot90(y_IN, 3);
    for j = 1:N+1
        stand_val = stand_val + coe(j) * x_sin .^ (j-1);
    end
end