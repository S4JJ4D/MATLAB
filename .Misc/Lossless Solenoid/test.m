%%
syms alpha beta real positive;

F = [7.5, .75];
x = 1e-3*[1, 5];
i0 = 0.05;

eq = [...
    2*F(1)/(i0^2) == beta/((alpha + beta*x(1))^2), ...
    2*F(2)/(i0^2) == beta/((alpha + beta*x(2))^2)];

res = vpasolve(eq);

alpha1 = double(res.alpha);
beta1 = double(res.beta);

f = 1/2 * i0^2 * beta1 ./ ((alpha1 + beta1 * x).^2)

f = 1/2 * i0^2 * beta0 ./ ((alpha0 + beta0 * x).^2)