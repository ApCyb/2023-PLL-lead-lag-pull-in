import numpy as np
from matplotlib import pyplot as plt
from datetime import datetime

from K_vco_ht_function import K_vco_ht_function
from K_vco_pt_function import K_vco_pt_function
from must_be_greater_than_1_by_pi import must_be_greater_than_1_by_pi
from omega_p_function import omega_p_function


def plot_bifurcation_diagram_3parameters(k, tau_1, tau_2):
    """
    Plots value of normalized pull-in frequency
    for slope k > 1/pi and values tau_1 > 0, tau_2 >= 0

    Parameters:
        k (float): must be greater than 1/pi
        tau_1 (float): must be positive
        tau_2 (float): must be non-negative
    """
    # Input validation
    must_be_greater_than_1_by_pi(k)
    if tau_1 <= 0:
        raise ValueError("tau_1 must be positive")
    if tau_2 < 0:
        raise ValueError("tau_2 must be non-negative")

    # Define the range of K_vco values for modeling
    K_vcos = np.concatenate([
        np.arange(0.1, 10.1, 0.1),
        np.arange(11, 101, 1),
        np.arange(110, 1001, 10),
        np.arange(1100, 10001, 100),
        np.arange(11000, 100001, 1000)
    ])

    n = len(K_vcos)  # Number of points
    omega_p = np.zeros(n)  # Initialize pull-in frequency array

    # Compute omega_p for each K_vco value
    for i in range(n):
        K_vco = K_vcos[i]
        omega_p[i] = omega_p_function(k, K_vco, tau_1, tau_2)

    # Normalize pull-in frequency
    omega_p_normalized = omega_p / K_vcos

    plt.figure(figsize=(10, 6))
    plt.semilogx(K_vcos, omega_p_normalized, 'black', linewidth=1)
    plt.grid(True)

    # Compute and plot the critical point x_ht
    # Critical points K_vco_ht and K_vco_pt
    K_vco_ht = K_vco_ht_function(k, tau_1, tau_2)

    # Plot black cross for K_vco_ht
    plt.semilogx(K_vco_ht, 1, 'kx', markersize=10, markeredgewidth=1)

    if tau_2 != 0:
        # Additional thresholds when tau_2 is non-zero
        K_vco_pt = K_vco_pt_function(k, tau_1, tau_2)
        omega_p_switching_normalized = omega_p_function(k, K_vco_pt, tau_1, tau_2)/K_vco_pt
        # Plot red cross for K_vco_pt
        plt.semilogx(K_vco_pt, omega_p_switching_normalized, 'rx', markersize=10, markeredgewidth=1)

    # Formatting the plot
    plt.xlabel(r'$K_{\rm vco}$', fontsize=20)
    plt.ylabel(r'$\frac{\omega_p}{K_{\rm vco}}$', fontsize=20)
    plt.xlim(min(K_vcos), max(K_vcos))
    plt.ylim(0, 1.1)

    plt.xticks([0.01, 0.1, 1, 10, 100, 1000, 5000],
               [r'$10^{-2}$', r'$10^{-1}$', r'$10^0$', r'$10^1$', r'$10^2$', r'$10^3$', r'$5\cdot10^3$'])

    plt.tick_params(axis='both', which='major', labelsize=15)

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    filename = f"{"bifurcation_3params"}_{timestamp}.pdf"
    plt.savefig(filename, format='pdf', bbox_inches='tight')
    plt.show()