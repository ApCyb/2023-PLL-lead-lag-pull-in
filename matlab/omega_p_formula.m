function omega_p = omega_p_formula(k, K_vco, tau_1, tau_2)
%OMEGA_P_FORMULA computes the pull-in frequency (omega_p) based on the system's parameters.
%
% Parameters:
%   k       - Slope of the system (dimensionless)
%   K_vco   - VCO (Voltage-Controlled Oscillator) gain
%   tau_1   - Time constant tau_1
%   tau_2   - Time constant tau_2
%
% Returns:
%   omega_p - Pull-in frequency based on different conditions

% Step 1: Compute the phase difference parameter and threshold values
mu = pi * k - 1;  
K_vco_ht = 1 / (k * (2 * tau_1 + tau_2 + 2 * sqrt(tau_1 * (tau_1 + tau_2)))); % Threshold for K_vco when tau_2 = 0

if (tau_2 == 0)
    % Case 1: When tau_2 is zero
    if (K_vco <= K_vco_ht)
        % Sub-case 1a: K_vco is less than or equal to the threshold
        omega_p = K_vco;
    else
        % Sub-case 1b: K_vco exceeds the threshold
        omega_p = omega_sep_formula(k, K_vco, tau_1, tau_2); % Use the separatrix formula
    end
else
    % Case 2: When tau_2 is non-zero
    k_ht = 2 * (tau_1 + tau_2 + sqrt(tau_1 * (tau_1 + tau_2))) / ...
           (pi * (2 * tau_1 + tau_2 + 2 * sqrt(tau_1 * (tau_1 + tau_2)))); % Threshold for k
    K_vco_pt = max(K_vco_ht, mu / (k * tau_2)); % Maximum allowable K_vco for the separatrix region

    if (k <= k_ht)
        % Sub-case 2a: k is less than or equal to the threshold k_ht
        if (K_vco <= K_vco_ht)
            % Sub-case 2a1: K_vco is less than or equal to the threshold
            omega_p = K_vco;
        else
            % Sub-case 2a2: K_vco exceeds the threshold
            omega_p = omega_ss_formula(k, K_vco, tau_1, tau_2); % Use the steady-state formula
        end
    else
        % Sub-case 2b: k exceeds the threshold k_ht
        if (K_vco <= K_vco_ht)
            % Sub-case 2b1: K_vco is less than or equal to the threshold
            omega_p = K_vco;
        else
            if (K_vco <= K_vco_pt)
                % Sub-case 2b2: K_vco is within the separatrix range
                omega_p = omega_sep_formula(k, K_vco, tau_1, tau_2); % Use the separatrix formula
            else
                % Sub-case 2b3: K_vco exceeds the separatrix range
                omega_p = omega_ss_formula(k, K_vco, tau_1, tau_2); % Use the steady-state formula
            end
        end
    end
end
end
