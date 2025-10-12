tol = 1e-10;
fx = @(x) atan(x) ./ x;
I = 1/5 * (-integral(fx, 0, 5, 'AbsTol', tol))

[n, err] = findN(tol, I)

function y = f(x)
    y = -1/5 * atan(x) ./ x;
end

function [n, err] = findN(tol, I)
    n = 2;
    while true
        syms x;
        P = legendreP(n,x);      
        x_nodes = sort(double(solve(P == 0)));

        a = zeros(size(x_nodes));
        for i = 1:length(x_nodes)
            others = x_nodes;
            others(i) = [];                       
            denom = prod(x_nodes(i) - others);     
            coeffs = poly(others);                

            p = numel(coeffs)-1:-1:0;
            integral_terms = coeffs ./ (p+1) .* (1 - (-1).^(p+1));
            a(i) = sum(integral_terms) / denom;
        end

        fvals = f(5/2 .* x_nodes + 5/2);
        integral_approx =  sum(a(:) .* fvals(:) *5/2);

        if abs(I - integral_approx) < tol
            err = abs(I - integral_approx);
            break;
        else
            n = n + 1;
        end

    end
end
