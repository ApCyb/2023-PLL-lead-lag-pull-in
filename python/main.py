import math

from plot_bifurcation_diagram_2parameters import plot_bifurcation_diagram_2parameters
from plot_bifurcation_diagram_3parameters import plot_bifurcation_diagram_3parameters

plot_bifurcation_diagram_2parameters(2/math.pi, [0, 0.005, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.7, 0.9])
plot_bifurcation_diagram_3parameters(k=0.5, tau_1=1.0, tau_2=0.5)
