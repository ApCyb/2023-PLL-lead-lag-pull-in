function K_vco_ht = K_vco_ht_function(k, tau_1, tau_2)
%K_vco_ht_function calculates K_vco_ht
 
K_vco_ht = 1 / (k*(2*tau_1 + tau_2 + 2*sqrt(tau_1*(tau_1 + tau_2))));
  
end

