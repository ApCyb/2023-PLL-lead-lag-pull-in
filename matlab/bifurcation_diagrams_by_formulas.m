% Close all figures, clear variables, and clear the command window
close all;
clc;
clear all;

% Define slope parameter k and parameter mu
k = 2 / pi;
mu = pi * k - 1;

% Array of a values (less than 1)
a_array = [0, 0.005, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.7, 0.9];
% Note: Standard engineering parameters tau_1 = 0.0448, tau_2 = 0.0185
% correspond to a = 0.2923.

% Array of x values for plotting
x_array = [0.0001:0.01:10, 11:0.1:100, 110:1:1000, 1100:10:10000];

% Loop over each value of a in the array
for i = 1:length(a_array)
    a = a_array(i); % Current value of a

    % Initialize array for normalized pull-in frequency y_p
    y_p = zeros(1, length(x_array));

    % Loop over each x value to compute y_p
    for j = 1:length(x_array)
        x = x_array(j);

        % Transformations:
        %   tau_1 = 1
        %   a = tau_2 / (1 + tau_2) => tau_2 = a / (1 - a)
        %   K_vco = x * (1 - a)
        tau_1 = 1;
        tau_2 = a / (1 - a);
        K_vco = x * (1 - a);

        % Compute normalized pull-in frequency
        y_p(j) = omega_p_formula(k, K_vco, tau_1, tau_2);
    end

    % Plot normalized pull-in frequency as a function of x
    semilogx(x_array, y_p ./ ((1 - a) * x_array), 'black', 'LineWidth', 1);
    grid on;
    hold on;

    % Compute and plot the critical point x_ht
    x_ht = 1 / (k * (2 - a + 2 * sqrt(1 - a))); % Threshold x_ht
    semilogx(x_ht, 1, 'x', 'LineWidth', 2, 'Color', 'black');

    % If a is nonzero, compute and plot x_pt
    if a ~= 0
        x_pt = max(mu / (a * k), x_ht); % Threshold x_pt
        y_pt = omega_p_formula(k, (1 - a) * x_pt, 1, a / (1 - a)) / ((1 - a) * x_pt);
        semilogx(x_pt, y_pt, 'x', 'LineWidth', 2, 'Color', 'red');
    end
end

% Set font size for the axes
set(gca, 'FontSize', 15);

% Label the axes with LaTeX formatting
xlabel('\textbf{$(\tau_1 + \tau_2)K_{\rm vco}$}', 'Interpreter', 'latex', 'fontsize', 20);
ylabel('\textbf{$\frac{\omega_l}{K_{\rm vco}}$}', 'Interpreter', 'latex', 'fontsize', 20);

% Adjust axis limits and ticks
axis([0.4, max(x_array), 0, 1.01]);
xticks([0.01, 0.1, 1, 10, 100, 1000, 5000]);
xticklabels({'10^{-2}', '10^{-1}', '10^0', '10^1', '10^2', '10^3', '5\cdot10^3'});
