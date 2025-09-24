function K_vco_pt = K_vco_pt_function(k, tau_1, tau_2)
%K_vco_ht_function calculates K_vco_ht
K_vco_ht = K_vco_ht_function(k, tau_1, tau_2);
K_vco_pt = max(K_vco_ht, (pi * k - 1) / (k * tau_2));
end

