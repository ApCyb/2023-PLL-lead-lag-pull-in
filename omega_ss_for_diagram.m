function omega_ss_for_diagram = omega_ss_for_diagram(a, x_array, k)

omega_ss_for_diagram = zeros(1, length(x_array));

for j=1:length(x_array)
    x = x_array(j);
    mu = pi*k - 1;
    eta = (a*k*x - mu)/(2*sqrt(x));
    kappa = sqrt(eta^2 + mu*k);

f = @(z1) (curve2(a, x, k, z0(a, x, k, z1), z1) - curve1(a, x, k, z0(a, x, k, z1), z1));

z1_start = eta + kappa;
z1_finish = k*sqrt(a*x);

if(abs(x -  mu/(a*k)) < 1/a) %заглушка
    omega_ss_for_diagram(j) = omega_sep_for_diagram(a, x, k);
else
    z1_pt = fzero(f, [z1_start+0.0000000000001, z1_finish]);
    s2 = (z0(a, x, k, z1_pt) + eta - kappa).*(z0(a, x, k, z1_pt) + eta + kappa)./...
    (z1_pt - eta + kappa)./(z1_pt - eta - kappa).*...
    ((z0(a, x, k, z1_pt) + eta + kappa).*(z1_pt - eta + kappa)./...
    (z1_pt - eta - kappa)./(z0(a, x, k, z1_pt) + eta - kappa)).^(eta/kappa);

omega_ss_for_diagram(j) = (sqrt(s2) - 1)./(sqrt(s2) + 1);
end
end
end












function z0 = z0(a, x, k, z1)
mu = pi*k - 1;
eta = (a*k*x - mu)/(2*sqrt(x));
xi = (a*k*x + 1)/(2*sqrt(x));
z0 = ((1 + mu)*k*z1 - 2*k*(mu*xi + eta))./((1 + mu)*k - 2*(xi - eta)*z1);
end


function curve1 = curve1(a, x, k, z0, z1)
mu = pi*k - 1;
eta = (a*k*x - mu)/(2*sqrt(x));
kappa = sqrt(eta^2 + mu*k);
curve1 = ((z0 + eta).^2 - kappa^2)./((z1 - eta).^2 - kappa^2).*...
    ((z0 + eta + kappa).*(z1 + kappa - eta)./(z0 + eta - kappa)./(z1 - eta - kappa)).^(eta/kappa);
end


function curve2 = curve2(a, x, k, z0, z1)
    xi = (a*k*x + 1)/(2*sqrt(x));
    rho = sqrt(abs(xi^2 - k));
    x_plus = 1/(k*(2 - a - 2*sqrt(1 - a)));
    if(x < x_plus)
        curve2 = (z0.^2 + 2*xi*z0 + k)./(z1.^2 - 2*xi*z1 + k) .*...
        exp(2*xi/rho*(atan(rho./(z0 + xi)) - atan((z1 - xi)/rho) + pi/2));
    end
    if(x == x_plus)
        curve2 = ((z0 + sqrt(k))./(z1 - sqrt(k)).*...
        exp(sqrt(k)/(z0 + sqrt(k)) + sqrt(k)/(z1 - sqrt(k)))).^2;
    end
    if(x > x_plus)
        curve2 = (z0 + xi - rho).*(z0 + xi + rho)./(z1 - xi + rho)./(z1 - xi - rho).*...
    ((z0 + xi + rho).*(z1 + rho - xi)./(z0 + xi - rho)./(z1 - xi - rho)).^(xi/rho);
    end
end




