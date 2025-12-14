N = 20;
[D,x] = cheb(N); x = x(2:N); % spectral differentiation matrix
nu = 2e-4;
L = D^2; L = nu*L(2:N,2:N); % 2nd-order differentiation
w = .53*x + .47*sin(-1.5*pi*x) - x; % use w = u-x to make BCs homogeneous
u = [1;w+x;-1];

h = 1/100; % time step

A  = h*L;                 
E  = expm(A);             
E2 = expm(A/2); 

w0 = w;

%-----------by contour integral---------
[f1, f2, f3, Q] = coe_counter(L, N, h);
[tt_c, uu_c] = ETDRK4(f1, f2, f3, Q, E, E2, u, h, w0, x);

%----------by directly compute-----------
[f1, f2, f3, Q] = coe_directly(L, N, h);
[tt_d, uu_d] = ETDRK4(f1, f2, f3, Q, E, E2, u, h, w0, x);

% --------------plot figure-----------------:
figure;
surf([1;x;-1],tt_c,uu_c'), lighting phong, axis tight
view([-45 60]), colormap(cool), light('col',[1 1 0],'pos',[-10 0 10])

figure;
surf([1;x;-1],tt_d,uu_d'), lighting phong, axis tight
view([-45 60]), colormap(cool), light('col',[1 1 0],'pos',[-10 0 10])

%----------------function----------------

%-------------RK coefficients()integral------------
function [f1, f2, f3, Q] = coe_counter(L, N, h)
    M = 32; % no. of points for resolvent integral
    r = 15*exp(1i*pi*((1:M)-.5)/M); % points along complex circle
    A = h*L;

    I = eye(N-1); Z = zeros(N-1);
    f1 = Z; f2 = Z; f3 = Z; Q = Z;
    for j = 1:M
        z = r(j);
        zIA = (z*I - A) \ I;
        
        Q     = Q     + h*zIA*(exp(z/2)-1);
        f1 = f1 + h*zIA*(-4-z+exp(z)*(4-3*z+z^2))/z^2;
        f2  = f2  + h*zIA*(2+z+exp(z)*(z-2))/z^2;
        f3 = f3 + h*zIA*(-4-3*z-z^2+exp(z)*(4-z))/z^2;
    end
    f1 = real(f1/M); f2 = real(f2/M); f3 = real(f3/M); Q = real(Q/M);
end

%-------------RK coefficients(directly)------------
function [f1, f2, f3, Q, E, E2] = coe_directly(L, N, h)
    I  = speye(N-1);

    A  = h*L;         
    A2 = A*A;
    A3 = A2*A;

    E  = expm(A);             
    E2 = expm(A/2);           

    Q = h * (A \ (E2 - I));
    f1 = h * (A3 \ ( -4*I - A + E*(4*I - 3*A + A2) ));
    f2  = h * (A3 \ (  2*I + A + E*(A - 2*I)       ));
    f3 = h * (A3 \ ( -4*I - 3*A - A2 + E*(4*I - A) ));
end

%--------------time step------------------
function [tt, uu] = ETDRK4(f1, f2, f3, Q, E, E2, u, h, w, x)
    uu = u; tt = 0;
    tmax = 70; nmax = round(tmax/h); nplt = floor((tmax/70)/h);
    for n = 1:nmax
        t = n*h;
        Nu = (w+x) - (w+x).^3;
        a = E2*w + Q*Nu;
        Na = a + x - (a+x).^3;
        b = E2*w + Q*Na;
        Nb = b + x - (b+x).^3;
        c = E2*a + Q*(2*Nb-Nu);
        Nc = c + x - (c+x).^3;
        w = E*w + f1*Nu + 2*f2*(Na+Nb) + f3*Nc;

        %------------store data for plot--------------
        if mod(n,nplt)==0
            u = [1;w+x;-1];
            uu = [uu,u]; tt = [tt,t];
        end
    end
end