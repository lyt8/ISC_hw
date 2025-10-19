t_end = 10;
y0 = 0.2;
x_plot = linspace(0, t_end, 1001);
y_plot = g(x_plot, y0);

figure; hold on; grid on
plot(x_plot, y_plot, 'Color','w', LineWidth=3); hold on;
xlabel('x'); ylabel('y');

dt = [0.1,0.5,1.0,1.5,2.0, 3.0];
for h = dt
    [fu_un, step_1, un] = forward_euler(y0, @f, h, t_end);
    plot(step_1, un, 'LineStyle', '--', 'Marker', 'o');
end

legend('y', 'h = 0.1', 'h = 0.5', 'h = 1.0', 'h = 1.5', 'h = 2.0', 'h = 3.0');


function y = g(x, y0)
    y = y0 .* exp(x) ./ (1 + y0 .* exp(x));
end

function dydt = f(y)
    dydt = y * (1 - y);
end

function [fu_un, step, un] = forward_euler(y0, f, h, t_end)
    un = zeros(ceil(t_end/h), 1);
    step = zeros(ceil(t_end/h), 1);
    un(1) = y0 / (1 + y0);
    t = 0;
    u_b = y0;
    i = 2;
    while t < t_end
        u_a = u_b + h * f(u_b);
        u_b = u_a;
        t = t + h; 

        un(i) = u_a;
        step(i) = t;
        i = i + 1;
    end
    fu_un = u_b; 
end

