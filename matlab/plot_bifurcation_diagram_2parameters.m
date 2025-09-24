function plot_bifurcation_diagram_2parameters(k, a_array)
%bifurcation_diagram_function plots bifurcation diagram for slope k > 1/pi
% and value a = tau_2/(tau_1 + tau_2) < 1 or array of values a_array
arguments
        k {mustBeGreaterThan1_by_Pi(k)}
        a_array {mustBeGreaterThanOrEqual(a_array, 0), mustBeLessThan(a_array, 1)}
end

% Array of x values for plotting
x_array = [0.0001:0.01:10, 11:0.01:100, 110:1:1000, 1100:10:10000];

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

        % Compute pull-in frequency
        y_p(j) = omega_p_function(k, K_vco, tau_1, tau_2);
    end


    % Normalize pull-in frequency
    y_p_normalized = y_p ./ ((1 - a) * x_array);

    % Plot normalized pull-in frequency as a function of x
    semilogx(x_array, y_p_normalized, 'black', 'LineWidth', 1);
    grid on;
    hold on;

    % Compute and plot the critical point x_ht
    x_ht = 1 / (k * (2 - a + 2 * sqrt(1 - a))); % Threshold x_ht
    semilogx(x_ht, 1, 'x', 'LineWidth', 2, 'Color', 'black');

    % If a is nonzero, compute and plot x_pt
    if a ~= 0
        x_pt = max((pi * k - 1) / (a * k), x_ht); % Threshold x_pt
        y_pt = omega_p_function(k, (1 - a) * x_pt, 1, a / (1 - a)) / ((1 - a) * x_pt);
        semilogx(x_pt, y_pt, 'x', 'LineWidth', 2, 'Color', 'red');
    end
end

% Set font size for the axes
set(gca, 'FontSize', 15);

% Label the axes with LaTeX formatting
xlabel('\textbf{$(\tau_1 + \tau_2)K_{\rm vco}$}', 'Interpreter', 'latex', 'fontsize', 20);
ylabel('\textbf{$\frac{\omega_p}{K_{\rm vco}}$}', 'Interpreter', 'latex', 'fontsize', 20);

% Adjust axis limits and ticks
axis([0.4, max(x_array), 0, 1.01]);
xticks([0.01, 0.1, 1, 10, 100, 1000, 10000]);
xticklabels({'10^{-2}', '10^{-1}', '10^0', '10^1', '10^2', '10^3', '10^4'});
