import math
from scipy.optimize import brentq

from omega_ht_function import omega_ht_function


def omega_pt_function(k, K_vco, tau_1, tau_2, mu, xi, eta, rho, kappa):
    """
    OMEGA_SS_FORMULA calculates frequency for the semistable cycle

    Parameters:
      k > 1/pi   - Slope
      K_vco > K_vco_pt  - VCO (Voltage-Controlled Oscillator) gain
      tau_1 > 0  - Time constant tau_1
      tau_2 > 0  - Time constant tau_2

    Returns:
      omega_pt - frequency corresponding to the birth of semistable cycle
    """

    # Step 1: Determine the range for z1
    z1_start = kappa + eta
    z1_finish = k * math.sqrt(tau_2 * K_vco)

    # Step 2: Find z1 using a numerical method
    z1_pt = find_z1(k, K_vco, tau_1, tau_2, mu, xi, eta, rho, kappa, z1_start, z1_finish)

    # Step 3: Calculate omega_pt based on the value of z1_pt
    if math.isclose(z1_pt, z1_start):
        # Case 1: z1_pt equals the starting point
        omega_pt = omega_ht_function(k, K_vco, tau_1, tau_2, mu, xi, eta, rho, kappa)  # Use separatrix formula
    else:
        # Case 2: General case where z1_pt does not equal the starting point
        z0_pt = z0_function(z1_pt, k, mu, xi, eta)    # Calculate corresponding z0
        numerator = (z0_pt + eta)**2 - kappa**2
        denominator = (z1_pt - eta)**2 - kappa**2
        term1 = numerator / denominator

        term2_1 = (z0_pt + kappa + eta)*(z1_pt + kappa - eta)/(z0_pt - kappa + eta)/(z1_pt - kappa - eta)
        term2 = (term2_1) ** (eta / kappa)
        s_pt = term1 * term2
        omega_pt = (math.sqrt(s_pt) - 1) / (math.sqrt(s_pt) + 1) * K_vco

    return omega_pt

def find_z1(k, K_vco, tau_1, tau_2, mu, xi, eta, rho, kappa, z1_start, z1_finish):
    """FIND_Z1 finds the value of z1 using a recursive or numerical method.

    Parameters:
        k, K_vco, tau_1, tau_2 - System parameters
        z1_start, z1_finish    - Search range for z1

    Returns:
        found_z1 - Value of z1 that satisfies the main_curve condition
    """

    if z1_finish - z1_start < 1e-12:  # since main_curve_z(z1_start) = +inf,
        # we cannot use root finding function at this point
        found_z1 = z1_start  # Return the starting value if range is negligible
    else:
        z1_middle = z1_start + (z1_finish - z1_start)/2  # Midpoint of the range
        middle_val = main_curve(z1_middle, k, K_vco, tau_1, tau_2, mu, xi, eta, rho, kappa)

        if middle_val < 0:
            # If main_curve is negative, recurse on the lower half
            found_z1 = find_z1(k, K_vco, tau_1, tau_2, mu, xi, eta, rho, kappa,
                              z1_start, z1_middle)
        else:
            # Otherwise, find the root using brentq
            main_curve_z = lambda z1: main_curve(z1, k, K_vco, tau_1, tau_2,
                                               mu, xi, eta, rho, kappa)

            try:
                found_z1 = brentq(main_curve_z, z1_middle, z1_finish,
                                 rtol=1e-12, xtol=1e-12)
            except ValueError:
                # Fallback if brentq fails
                found_z1 = z1_middle

    return found_z1


def main_curve(z1, k, K_vco, tau_1, tau_2, mu, xi, eta, rho, kappa):
    """MAIN_CURVE calculates the difference between the left-hand and right-hand sides of the main equation.

    Parameters:
        z1 - Variable for the main curve equation
        k, K_vco, tau_1, tau_2 - System parameters

    Returns:
        main_curve_val - Value of the main curve function
    """

    z0 = z0_function(z1, k, mu, xi, eta)

    # Define thresholds for k and K_vco
    k_fn = 2*(tau_1 + tau_2 - math.sqrt(tau_1*(tau_1 + tau_2))) / (math.pi*(2*tau_1 + tau_2 - 2*math.sqrt(tau_1*(tau_1 + tau_2))))
    K_vco_fn = 1/(k*(2*tau_1 + tau_2 - 2*math.sqrt(tau_1*(tau_1 + tau_2))))

    # Calculate the right-hand side of the equation
    if k < k_fn and K_vco < K_vco_fn:
        numerator = z0**2 + 2*xi*z0 + k
        denominator = z1**2 - 2*xi*z1 + k
        atan_term = math.atan(rho/(z0 + xi)) - math.atan((z1 - xi)/rho) + math.pi/2
        right_hand_side = (numerator / denominator) * math.exp(2*xi/rho * atan_term)
    elif k < k_fn and math.isclose(K_vco, K_vco_fn):
        numerator = z0 + math.sqrt(k)
        denominator = z1 - math.sqrt(k)
        exp_term = math.exp(math.sqrt(k)/(z0 + math.sqrt(k)) + math.sqrt(k)/(z1 - math.sqrt(k)))
        right_hand_side = (numerator / denominator * exp_term)**2
    elif (k < k_fn and K_vco > K_vco_fn) or (k >= k_fn):
        term1 = (z0 + xi - rho)*(z0 + xi + rho) / ((z1 - xi + rho)*(z1 - xi - rho))
        term2_numerator = (z0 + xi + rho)*(z1 + rho - xi)
        term2_denominator = (z0 + xi - rho)*(z1 - xi - rho)
        term2 = (term2_numerator / term2_denominator) ** (xi/rho)
        right_hand_side = term1 * term2

    # Calculate the left-hand side of the equation
    numerator = (z0 + eta)**2 - kappa**2
    denominator = (z1 - eta)**2 - kappa**2
    term1 = numerator / denominator

    term2_numerator = (z0 + eta + kappa)*(z1 + kappa - eta)
    term2_denominator = (z0 + eta - kappa)*(z1 - eta - kappa)
    term2 = (term2_numerator / term2_denominator) ** (eta/kappa)

    left_hand_side = term1 * term2

    # Compute the main curve value as the difference
    main_curve_val = left_hand_side - right_hand_side

    return main_curve_val


def z0_function(z1, k, mu, xi, eta):
    """Z0 calculates the value of z0 based on z1 and system parameters.

    Parameters:
        z1, k, mu, xi, eta - Input parameters

    Returns:
        z0_val - Corresponding z0 value
    """

    numerator = (1 + mu)*k*z1 - 2*(mu*xi + eta)*k
    denominator = (1 + mu)*k - 2*(xi - eta)*z1
    z0_val = numerator / denominator

    return z0_val