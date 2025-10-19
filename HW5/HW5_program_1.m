x_plot = linspace(0, 10, 1001);
y_plot = g(x_plot);
dt = [0.4, 0.41, 0.1];
y0 = 1;

%----------forward euler method------------
figure; hold on; grid on
plot(x_plot, y_plot);
xlabel('t'); ylabel('y');
for h = dt
    [fu_un, step_1, un] = forward_euler(y0, @f, h);
    plot(step_1, un, 'LineStyle', '--', 'Marker', 'o');
end

legend('y', 'h = 0.4', 'h = 0.41', 'h = 0.1');
title('Forward Euler');

%------------back euler method--------------

for h = dt
    figure; hold on; grid on
    plot(x_plot, y_plot);
    xlabel('t'); ylabel('y');
    [bu_un, step_1, un] = back_euler(y0, h);
    plot(step_1, un, 'LineStyle', '--', 'Marker', 'o');
    legend('y', sprintf('h = %.2f', h))
    title(sprintf('Forward Euler with h = %.2f', h));
end

function y = g(x)
    y = exp(-5*x);
end

function dydt = f(y)
    dydt = -5 * y;
end

function [fu_un, step, un] = forward_euler(y0, f, h)
    un = zeros(ceil(10/h), 1);
    step = zeros(ceil(10/h), 1);
    un(1) = y0;
    t = 0;
    u_b = y0;
    i = 2;
    while t < 10
        u_a = u_b + h * f(u_b);
        u_b = u_a; 
        t = t + h;

        un(i) = u_b; 
        step(i) = t; 
        i = i + 1;
    end
    fu_un = u_b; 
end

function [bu_en, step, un] = back_euler(y0, h)
    un = zeros(ceil(10/h), 1);
    step = zeros(ceil(10/h), 1);
    un(1) = y0;
    t = 0;
    u_b = y0;
    i = 2;
    while t < 10
        u_a = u_b / (h*(-5)+1);
        u_b = u_a; 
        t = t + h; 

        un(i) = u_b; 
        step(i) = t; 
        i = i + 1;
    end
    bu_en = u_b; 
end