x_plot = linspace(-1.5, 1.5);
mmax = 20;
tol = 1e-13;
R = mean(f(x_plot));

Z = x_plot(:);
M = length(Z);  
F = f(Z);                
C = [];                   
J = 1:M;          
SF = spdiags(F,0, M,M);   
support_point = [];       
data_value = [];          
errvec = [];

for m = 1:mmax
    [~, index] = max(abs(F - R));
    new_point = Z(index);
    f_new = F(index);
    support_point = [support_point, new_point];
    data_value = [data_value, f_new];
    J(J == index) = [];
    C = [C, 1./(Z - new_point)];
    Sf = diag(data_value);
    A = SF*C - C*Sf;
    [~,~,V] = svd(A(J,:), 0);
    w = V(:, end);
    N = C*(w .* data_value.');
    D = C*w;
    R = F;
    R(J) = N(J)./D(J);

    err = norm(F - R, inf);
    errvec = [errvec; err];
    fprintf('m = %d, err = %.3e\n', m, err);

    if err <= tol * norm(F, inf)
        break;
    end
end

r = @(zz) get_rational_approx(zz, support_point, data_value, w); 
[pol,res] = prz(r,support_point, data_value, w);


xx = linspace(-3.5, 4.5,1005);
figure; 
plot(x_plot, f(x_plot), 'r', 'LineWidth',4); hold on;
plot(x_plot, r(x_plot), 'LineStyle', '--', 'Color', 'b', 'LineWidth',4); hold on;
title('approximate to gamma by AAA in [-1.5, 1.5]');
plot(support_point,r(support_point),'o','Color','k');
legend('gamma', 'AAA');

set(gcf, 'Color', 'w');
set(gca, 'Color', 'w');
set(gca, 'GridColor', [0.3 0.3 0.3]);  
set(gca, 'XColor', 'k', 'YColor', 'k', 'FontSize', 12, 'LineWidth', 1);
set(findall(gca, 'Type', 'text'), 'Color', 'k');  

figure;
plot(xx, f(xx), 'LineWidth',4); hold on;
plot(xx, r(xx), 'LineStyle', '--', 'LineWidth',4); hold on;
title('Extrapolated approximation over [-3.5, 4.5]');
legend('gamma', 'AAA');

set(gcf, 'Color', 'w');
set(gca, 'Color', 'w');
set(gca, 'GridColor', [0.3 0.3 0.3]);  
set(gca, 'XColor', 'k', 'YColor', 'k', 'FontSize', 12, 'LineWidth', 1);
set(findall(gca, 'Type', 'text'), 'Color', 'k');   

%--------------------function--------------
function y = f(x_plot)
    y = gamma(x_plot);
end

function r = get_rational_approx(zz,z,f,w)
    zv = zz(:);
    CC = 1 ./ (zv - z);
    r = (CC*(w.*f(:)))./(CC*w);
    ii = find(isnan(r));
    for j = 1:length(ii)
        r(ii(j)) = f(find(zv(ii(j))==z));
    end
    r = reshape(r, size(zz));
end

function [pol,res] = prz(r,z,f,w)
    m = length(w); 
    B = eye(m+1); B(1,1) = 0;
    
    % Poles
    E = [0 w.'; ones(m,1) diag(z)];
    pol = eig(E,B); 
    pol = pol(~isinf(pol));
    
    % Residues
    dz = 1e-5 * exp(2i*pi*(1:4)/4);
    res = r(bsxfun(@plus, pol, dz))*dz.'/4; 

end
    