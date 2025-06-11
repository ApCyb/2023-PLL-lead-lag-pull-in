import math
import numpy as np
from matplotlib import pyplot as plt
from datetime import datetime

from must_be_greater_than_1_by_pi import must_be_greater_than_1_by_pi
from omega_p_function import omega_p_function

def plot_bifurcation_diagram_2parameters(k, a_array):
    """
    Plots bifurcation diagram for slope k > 1/pi and value(s) a = tau_2/(tau_1 + tau_2) < 1

    Parameters:
        k (float): must be greater than 1/pi
        a_array (array-like): value(s) in [0, 1)
    """
    # Input validation
    must_be_greater_than_1_by_pi(k)
    a_array = np.asarray(a_array)
    if np.any(a_array < 0) or np.any(a_array >= 1):
        raise ValueError("a_array must be in [0, 1)")

    # Create x array with same ranges as MATLAB version
    x_ranges = [
        np.arange(0.0001, 10, 0.01),
        np.arange(11, 100, 0.01),
        np.arange(110, 1000, 1),
        np.arange(1100, 10000, 10)
    ]
    x_array = np.concatenate(x_ranges)

    # Setup figure
    plt.figure(figsize=(10, 6))
    plt.grid(True)

    # Loop through each a value
    for a in a_array:
        # Initialize array for pull-in frequency
        y_p = np.zeros(len(x_array))

        # Calculate pull-in frequency for each x value
        for j, x in enumerate(x_array):
            # Parameter transformations:
            tau_1 = 1.0
            tau_2 = a / (1 - a)
            K_vco = x * (1 - a)

            # Compute pull-in frequency using implemented function
            y_p[j] = omega_p_function(k, K_vco, tau_1, tau_2)

        # Normalize pull-in frequency
        y_p_normalized = y_p / ((1 - a) * x_array)

        # Plot normalized pull-in frequency vs x (log scale)
        plt.semilogx(x_array, y_p_normalized, 'black', linewidth=1)

        # Compute and plot heteroclinic threshold (x_ht)
        x_ht = 1 / (k * (2 - a + 2 * math.sqrt(1 - a)))
        plt.semilogx(x_ht, 1, 'kx', markersize=10, markeredgewidth=1)

        # If a ≠ 0, compute and plot pull-in threshold (x_pt)
        if a != 0:
            x_pt = max((math.pi * k - 1) / (a * k), x_ht)
            y_pt = omega_p_function(k, (1 - a) * x_pt, 1, a / (1 - a)) / ((1 - a) * x_pt)
            plt.semilogx(x_pt, y_pt, 'rx', markersize=10, markeredgewidth=1)
    # Format plot
    plt.xlabel(r'$(\tau_1 + \tau_2)K_{\rm vco}$', fontsize=20)
    plt.ylabel(r'$\frac{\omega_p}{K_{\rm vco}}$', fontsize=20)
    plt.xlim(0.4, 10000)
    plt.ylim(0, 1.01)
    plt.xticks([0.01, 0.1, 1, 10, 100, 1000, 10000],
               ['$10^{-2}$', '$10^{-1}$', '$10^0$', '$10^1$', '$10^2$', '$10^3$', '$10^4$'])
    plt.tick_params(axis='both', which='major', labelsize=15)

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    filename = f"{"bifurcation_2params"}_{timestamp}.pdf"
    plt.savefig(filename, format='pdf', bbox_inches='tight')
    plt.show()