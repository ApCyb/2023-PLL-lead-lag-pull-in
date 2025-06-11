import math
from K_vco_ht_function import K_vco_ht_function

def K_vco_pt_function(k, tau_1, tau_2):
    K_vco_ht = K_vco_ht_function(k, tau_1, tau_2)
    K_vco_pt_1 = (math.pi * k - 1) / (k * tau_2)
    return max(K_vco_ht, K_vco_pt_1)