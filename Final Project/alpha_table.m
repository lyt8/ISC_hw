h = 1;
R = 2.0;
M = 32;

zlist = [10, 1, 1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9];

f1 = @(z) alpha_fun(z, h);
f2 = @(z) beta_fun(z, h);

fprintf('%8s %22s %22s %22s\n', 'z', 'Formula', 'Contour(R=2,M=32)', 'Exact(vpa)');

for lam = zlist
    val_formula = f1(lam);                       
    val_contour = contour_avg_full(f1, lam, M, R); 
    val_exact   = exact_vpa(f1, lam, 120);       

    fprintf('%8.0e %22.14e %22.14e %22.14e\n', ...
        lam, val_formula, real(val_contour), val_exact);
end

fprintf('%8s %22s %22s %22s\n', 'z', 'Formula', 'Contour(R=2,M=32)', 'Exact(vpa)');

for lam = zlist
    val_formula = f2(lam);                       
    val_contour = contour_avg_full(f2, lam, M, R); 
    val_exact   = exact_vpa(f2, lam, 120);       

    fprintf('%8.0e %22.14e %22.14e %22.14e\n', ...
        lam, val_formula, real(val_contour), val_exact);
end
% ===================== COEFFICIENT DEFINITION =====================
function a = alpha_fun(z, h)
    a = h^-2 .* ( -4 - z*h + exp(z*h).*(4 - 3.*z*h + (z*h).^2) )/ (z.^3);
end

function b = beta_fun(z, h)
    b = h^-2 .*( 2 + z.* h + exp(z*h).*(-2 + z*h))/ (z.^3);
end

function val = contour_avg_full(f, lam, M, R)
    val = 0;
    theta = 2*pi*((1:M)-0.5)/M;   
    r = lam + R*exp(1i*theta);
    for j = 1:M
        z = r(j);
       val =val + f(z);
    end
    val = val / M;
end

function val = exact_vpa(f, lam, ndigits)
    if nargin < 3, ndigits = 80; end
    digits(ndigits);
    lams = vpa(lam);
    vals = f(lams);
    val  = double(vals);
end
