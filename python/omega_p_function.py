import math

from K_vco_ht_function import K_vco_ht_function
from K_vco_pt_function import K_vco_pt_function
from must_be_greater_than_1_by_pi import must_be_greater_than_1_by_pi
from omega_ht_function import omega_ht_function
from omega_pt_function import omega_pt_function


def omega_p_function(k, K_vco, tau_1, tau_2):
    """
    OMEGA_P_FORMULA computes the pull-in frequency (omega_p) based on the system's parameters.

    Parameters:
      k > 1/pi   - Slope
      K_vco > 0  - VCO (Voltage-Controlled Oscillator) gain
      tau_1 > 0  - Time constant tau_1
      tau_2 >= 0  - Time constant tau_2

    Returns:
      omega_p - Pull-in frequency based on different conditions
    """

    must_be_greater_than_1_by_pi(k)
    if K_vco <= 0:
        raise ValueError("K_vco must be positive.")
    if tau_1 <= 0:
        raise ValueError("tau_1 must be positive.")
    if tau_2 < 0:
        raise ValueError("tau_2 must be non-negative.")

    # Step 1: Compute the phase difference parameter and threshold values
    K_vco_ht = K_vco_ht_function(k, tau_1, tau_2)

    # Step 2: Compute intermediate parameters for omega_ht and omega_pt
    mu = math.pi * k - 1
    xi = (1 + k * tau_2 * K_vco) / (2 * math.sqrt((tau_1 + tau_2) * K_vco))
    eta = (k * tau_2 * K_vco - mu) / (2 * math.sqrt((tau_1 + tau_2) * K_vco))
    rho = math.sqrt(abs(xi**2 - k))
    kappa = math.sqrt(eta**2 + k * mu)

    if tau_2 == 0:
        # Case 1: When tau_2 is zero
        if K_vco <= K_vco_ht:
            # Sub-case 1a: K_vco is less than or equal to the threshold
            omega_p = K_vco
        else:
            # Sub-case 1b: K_vco exceeds the threshold
            omega_p = omega_ht_function(k, K_vco, tau_1, tau_2, mu, xi, eta, rho, kappa)
    else:
        # Case 2: When tau_2 is non-zero
        k_ht = 2 * (tau_1 + tau_2 + math.sqrt(tau_1 * (tau_1 + tau_2))) / (
               math.pi * (2 * tau_1 + tau_2 + 2 * math.sqrt(tau_1 * (tau_1 + tau_2))))
        K_vco_pt = K_vco_pt_function(k, tau_1, tau_2)

        if k <= k_ht:
            # Sub-case 2a: k is less than or equal to the threshold k_ht
            if K_vco <= K_vco_ht:
                # Sub-case 2a1: K_vco is less than or equal to the threshold
                omega_p = K_vco
            else:
                # Sub-case 2a2: K_vco exceeds the threshold
                omega_p = omega_pt_function(k, K_vco, tau_1, tau_2, mu, xi, eta, rho, kappa)
        else:
            # Sub-case 2b: k exceeds the threshold k_ht
            if K_vco <= K_vco_ht:
                # Sub-case 2b1: K_vco is less than or equal to the threshold
                omega_p = K_vco
            else:
                if K_vco <= K_vco_pt:
                    # Sub-case 2b2: K_vco is within the separatrix range
                    omega_p = omega_ht_function(k, K_vco, tau_1, tau_2, mu, xi, eta, rho, kappa)
                else:
                    # Sub-case 2b3: K_vco exceeds the separatrix range
                    omega_p = omega_pt_function(k, K_vco, tau_1, tau_2, mu, xi, eta, rho, kappa)

    return omega_p
