function plot_bifurcation_diagram_3parameters(k, tau_1, tau_2)
%omega_p_normalized_function plots value of normalized pull-in frequency
% for slope k > 1/pi and values tau_1 > 0, tau_2 >= 0
arguments
        k {mustBeGreaterThan1_by_Pi(k)}
        tau_1 {mustBePositive}
        tau_2 {mustBeNonnegative}
end

% Define the range of K_vco values for modeling
K_vcos = [0.1:0.1:10, 11:1:100, 110:10:1000, 1100:100:10000, 11000:1000:100000];
n = length(K_vcos); % Number of points
omega_p = zeros(1, n); % Initialize pull-in frequency array

% Compute omega_p for each K_vco value
for i = 1:n
    K_vco = K_vcos(i);
    omega_p(i) = omega_p_function(k, K_vco, tau_1, tau_2);
end

% Normalize pull-in frequency
omega_p_normalized = omega_p ./ K_vcos;

% Plot normalized pull-in frequency (omega_p / K_vco) on a semilogarithmic scale
semilogx(K_vcos, omega_p_normalized, 'black', 'LineWidth', 1);
 
grid on;
hold on;

 % Compute and plot the critical point x_ht
 % Critical points K_vco_ht and K_vco_pt
K_vco_ht = K_vco_ht_function(k, tau_1, tau_2);
if tau_2 ~= 0
    % Additional thresholds when tau_2 is non-zero
    K_vco_pt = K_vco_pt_function(k, tau_1, tau_2);
end

semilogx(K_vco_ht, 1, 'x', 'LineWidth', 2, 'Color', 'black');
omega_p_switching_normalized = omega_p_function(k, K_vco_pt, tau_1, tau_2)/K_vco_pt;
semilogx(K_vco_pt, omega_p_switching_normalized, 'x', 'LineWidth', 2, 'Color', 'red');

% Formatting the plot
set(gca, 'FontSize', 15);
xlabel('\textbf{$K_{\rm vco}$}', 'Interpreter', 'latex', 'fontsize', 20);
ylabel('\textbf{$\frac{\omega_p}{K_{\rm vco}}$}', 'Interpreter', 'latex', 'fontsize', 20);
axis([min(K_vcos), max(K_vcos), 0, 1.1]); % Set axis limits
xticks([0.01, 0.1, 1, 10, 100, 1000, 5000]);
xticklabels({'10^{-2}', '10^{-1}', '10^0', '10^1', '10^2', '10^3', '5\cdot10^3'});


 
