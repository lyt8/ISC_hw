x_co = linspace(-5, 5, 2001);
y_co = f(xx);

N = 15;
x = linspace(-5,5,N+1);
y = f(x);

w = barycentric_weights(x, N);
m = barycentric_value(x, y, x_co, w);

plot(xx, yy); hold on
plot(xx, m);
xlabel('x'); ylabel('y')
title('Runge f(x)=1/(1+x^2)')
legend('f(x)', 'Barycentric')

function y = f(x)
    y = 1 ./ (1 + x.^2);
end

function w = barycentric_weights(x, N)
w = ones(1, N+1);
for i = 1:N+1
    for j = [1:i-1, i+1:N+1]
        w(i) = w(i) / (x(i) - x(j));
    end
end
end

function m = barycentric_value(x, y, x_co, w)
m = ones(1, 2001);
for i = 1:2001
     [dmin, j] = min(abs(x_co(i) - x));
     if dmin < 1e-6
         m(i) = y(j);      
     else
         t = w ./ (x_co(i) - x);          
         m(i) = sum(t .* y) / sum(t); 
     end
end
end

