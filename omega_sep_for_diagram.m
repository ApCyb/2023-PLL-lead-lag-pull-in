function omega_sep_for_diagram = omega_sep_for_diagram(a, x, k)

mu = pi*k - 1;
xi = (a*k*x + 1)./(2*sqrt(x));
eta = (a*k*x - mu)./(2*sqrt(x));
rho = sqrt(abs(xi.^2 - k));
kappa = sqrt(eta.^2 + mu*k);



x_ht = 1/(k*(2 - a + 2*sqrt(1 - a)));
x_nf = 1/(k*(2 - a - 2*sqrt(1 - a)));

s = ((kappa - eta).^2 + 2*xi.*(kappa - eta) + k)./((kappa + eta).^2 - 2*xi.*(kappa + eta) + k).*...
            exp(2*xi./rho.*(atan(((xi - eta).^2 + rho.^2 - kappa.^2)./(2*rho.*kappa)) + pi/2)) .* (x > x_ht) .* (x <= x_nf) +...
((kappa - eta + xi).^2 - rho.^2)./((kappa + eta - xi).^2 - rho.^2).*...
        ((((kappa + rho).^2) - (xi - eta).^2)./((kappa - rho).^2 - (xi - eta).^2)).^(xi./rho) .* (x > x_nf);
    
    omega_sep_for_diagram = (sqrt(s) - 1)./(sqrt(s) + 1);
end


