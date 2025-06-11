import math

def omega_ht_function(k, K_vco, tau_1, tau_2, mu, xi, eta, rho, kappa):
    """
    omega_ht_function calculates frequency for the separatrix cycle

    Parameters:
      k > 1/pi   - Slope
      K_vco > K_vco_ht  - VCO (Voltage-Controlled Oscillator) gain
      tau_1 > 0  - Time constant tau_1
      tau_2 >= 0  - Time constant tau_2

    Returns:
      omega_ht - frequency corresponding to heteroclinic trajectory
    """

    # Threshold K_vco_fn is only defined if tau_2 is nonzero
    K_vco_fn = None
    if tau_2 != 0:
        K_vco_fn = 1 / (k * (2 * tau_1 + tau_2 - 2 * math.sqrt(tau_1 * (tau_1 + tau_2))))

    # Step 1: Calculate s_ht
    if tau_2 == 0 or (K_vco_fn is not None and K_vco < K_vco_fn):
        # Case 1: tau_2 == 0 or K_vco is below the threshold
        numerator = (kappa - eta)**2 + 2 * xi * (kappa - eta) + k
        denominator = (kappa + eta)**2 - 2 * xi * (kappa + eta) + k
        atan_arg = ((xi - eta)**2 + rho**2 - kappa**2) / (2 * rho * kappa)
        exp_term = math.exp((2 * xi / rho) * (math.atan(atan_arg) + math.pi / 2))
        s_ht = (numerator / denominator) * exp_term

    elif K_vco_fn is not None and math.isclose(K_vco, K_vco_fn):
        # Case 2: K_vco equals the threshold
        sqrt_k = math.sqrt(k)
        numerator = kappa - eta + sqrt_k
        denominator = ((kappa + eta - sqrt_k) * math.exp(2 * sqrt_k * kappa / (kappa**2 - (eta - sqrt_k)**2)))**2
        s_ht = (numerator / denominator)

    elif K_vco_fn is not None and K_vco > K_vco_fn:
        # Case 3: K_vco is above the threshold
        term1 = ((kappa - eta + xi)**2 - rho**2) / ((kappa + eta - xi)**2 - rho**2)
        term2_numerator = (kappa + rho)**2 - (xi - eta)**2
        term2_denominator = (kappa - rho)**2 - (xi - eta)**2
        term2 = (term2_numerator / term2_denominator) ** (xi / rho)
        s_ht = term1 * term2
    else:
        raise ValueError("Incorrect parameters")

    # Step 2: Calculate the frequency
    sqrt_s_ht = math.sqrt(s_ht)
    omega_ht = (sqrt_s_ht - 1) / (sqrt_s_ht + 1) * K_vco

    return omega_ht