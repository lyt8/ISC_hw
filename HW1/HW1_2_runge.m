x_ru = linspace(-5, 5, 2001);
y_ru = f(x_ru);

N = 15;
x_IN = linspace(-5,5,N+1);
y_IN = f(x_IN);

w = barycentric_weights(x_IN, N);
m = barycentric_value(x_IN, y_IN, x_ru, w);

plot(x_ru, y_ru); hold on
plot(x_ru, m);
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

function m = barycentric_value(x_IN, y_IN, x_ru, w)
m = ones(1, 2001);
for i = 1:2001
     [dmin, j] = min(abs(x_ru(i) - x_IN));
     if dmin < 1e-6
         m(i) = y_IN(j);      
     else
         t = w ./ (x_ru(i) - x_IN);          
         m(i) = sum(t .* y_IN) / sum(t); 
     end
end
end

