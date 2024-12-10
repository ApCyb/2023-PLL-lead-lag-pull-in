function omega_ht = omega_ht_function(k, K_vco, tau_1, tau_2)
%omega_ht_function calculates frequency for the separatrix cycle
%
% Parameters:
%   k       - slope
%   K_vco   - VCO gain
%   tau_1   - Time constant tau_1
%   tau_2   - Time constant tau_2
%
% Returns:
%   omega_ht - frequency corresponding to heteroclinic trakectory

% Step 1: Compute intermediate parameters
mu = pi*k - 1;
xi = (1 + k*tau_2*K_vco)/(2*sqrt((tau_1+tau_2)*K_vco));
eta = (k*tau_2*K_vco - mu)/(2*sqrt((tau_1+tau_2)*K_vco));
rho = sqrt(abs(xi^2 - k));
kappa = sqrt(eta^2 + k*mu);

% Step 2: Compute critical VCO gain thresholds
K_vco_ht = 1 / (k*(2*tau_1 + tau_2 + 2*sqrt(tau_1*(tau_1 + tau_2))));
% Threshold K_vco_fn is only defined if tau_2 is nonzero
if(tau_2 ~= 0)
    K_vco_fn = 1 / (k*(2*tau_1 + tau_2 - 2*sqrt(tau_1*(tau_1 + tau_2))));
end


% Step 3: Determine the separatrix cycle existence conditions and calculate s_ht
if(tau_2 ==0 || K_vco < K_vco_fn)
    % Case 1: tau_2 == 0 or K_vco is below the threshold
    s_ht = ((kappa - eta)^2 + 2*xi*(kappa - eta) + k)/((kappa + eta)^2 - 2*xi*(kappa + eta) + k)*...
            exp(2*xi/rho*(atan(((xi - eta)^2 + rho^2 - kappa^2)/(2*rho*kappa)) + pi/2));
else
    if(K_vco == K_vco_fn)
    % Case 2: K_vco equals the threshold
        s_ht = (kappa - eta + sqrt(k))/((kappa + eta - sqrt(k))*...
        exp(2*sqrt(k)*kappa/(kappa^2 - (eta - sqrt(k))^2)))^2;
    else
        if(K_vco > K_vco_fn)
        % Case 3: K_vco is above the threshold
            s_ht = ((kappa - eta + xi)^2 - rho^2)/((kappa + eta - xi)^2 - rho^2)*...
        ((((kappa + rho)^2) - (xi - eta)^2)/((kappa - rho)^2 - (xi - eta)^2))^(xi/rho);
        end
    end
end

% Step 4: Calculate the frequency
omega_ht = (sqrt(s_ht) - 1)/(sqrt(s_ht) + 1)*K_vco;
end

