import math

def K_vco_ht_function(k, tau_1, tau_2):
    K_vco_ht = 1 / (k * (2 * tau_1 + tau_2 + 2 * math.sqrt(tau_1 * (tau_1 + tau_2))))
    return K_vco_ht