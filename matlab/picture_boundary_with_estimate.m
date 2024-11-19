% Close all figures, clear variables, and clear command window
close all;
clc;

% Define the slope parameter k
k = 2 / pi;

% Engineering parameters
tau_1 = 0.0448; % Time constant tau_1
tau_2 = 0.0185; % Time constant tau_2

mu = pi * k - 1;

% Critical points for switching formulas
K_vco_ht = 1 / (k * (2 * tau_1 + tau_2 + 2 * sqrt(tau_1 * (tau_1 + tau_2)))); % Upper threshold
if tau_2 ~= 0
    % Additional thresholds when tau_2 is non-zero
    K_vco_fn = 1 / (k * (2 * tau_1 + tau_2 - 2 * sqrt(tau_1 * (tau_1 + tau_2))));
    K_vco_pt = max(mu / (k * tau_2), K_vco_ht); % Maximum threshold for separatrix
end

% Define the range of K_vco values for modeling
K_vcos = [0.1:0.1:10, 11:1:100, 110:10:1000, 1100:100:10000, 11000:1000:100000];
n = length(K_vcos); % Number of points
omega_p = zeros(1, n); % Initialize pull-in frequency array

% Compute omega_p for each K_vco value
for i = 1:n
    K_vco = K_vcos(i);
    omega_p(i) = omega_p_formula(k, K_vco, tau_1, tau_2);
end

% Plot normalized pull-in frequency (omega_p / K_vco) on a semilogarithmic scale
semilogx(K_vcos, omega_p ./ K_vcos, 'black', 'LineWidth', 3);
grid on;
hold on;

if tau_2 ~= 0
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Plot heteroclinic trajectory for K_vco > K_vco_pt
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Select K_vco values greater than the K_vco_pt threshold
    K_vcos_semistable_cycle = K_vcos(K_vcos > K_vco_pt);
    m = length(K_vcos_semistable_cycle); % Number of points
    omega_sep = zeros(1, m); % Initialize omega_sep array

    % Compute omega_sep for each K_vco
    for i = 1:m
        K_vco = K_vcos_semistable_cycle(i);
        omega_sep(i) = omega_sep_formula(k, K_vco, tau_1, tau_2);
    end
    % Plot separatrix trajectory
    semilogx(K_vcos_semistable_cycle, omega_sep ./ K_vcos_semistable_cycle, 'blue--', 'LineWidth', 1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute and plot the asymptotic formula
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a = tau_2 / (tau_1 + tau_2); % Parameter for asymptotics
if a ~= 0
    % Define the function to solve for b (root-finding)
    fcn = @(x) a * (2 * x - a - x^2) / (x * (x - a)) - log((x^2 * (1 - a) / (x - a)^2));
    start = a + eps; % Start point for root-finding
    finish = sqrt(a); % End point
    b = fzero(fcn, [start, finish]); % Find the root
    y_as = (-2 * a * b + b^2 + a) / (2 * b - b^2 - a); % Asymptotic value

    % Lyapunov function pull-in estimate
    omega_Lyapunov = (tau_1 / (2 * sqrt(tau_2 * (tau_1 + tau_2)) - 2 * tau_2) - ...
                      sqrt(tau_1^2 / (2 * sqrt(tau_2 * (tau_1 + tau_2)) - 2 * tau_2)^2 - 1));
    % Plot Lyapunov estimate
    semilogx(K_vcos, omega_Lyapunov*K_vcos./K_vcos, 'red', 'LineWidth', 1);
else
    y_as = 0; % Set y_as to zero if a is zero
end
% Plot the horizontal line for the asymptotic value
semilogx([min(K_vcos), max(K_vcos)], [y_as, y_as], 'g', 'LineWidth', 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot vertical lines for critical K_vco values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
y = 0:0.0001:2; % Range for vertical lines
semilogx(K_vco_ht * y ./ y, y, 'black--', 'LineWidth', 1); % Vertical line for K_vco_ht
if tau_2 ~= 0
    semilogx(K_vco_pt * y ./ y, y, 'black--', 'LineWidth', 2); % Vertical line for K_vco_pt
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Lock-in %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% omega_l = zeros(1, n);
% for i = 1:n
%     K_vco = K_vcos(i);
%     omega_l(i) = omega_l_formula(tau_1, tau_2, k, K_vco);
% end
% loglog(1./K_vcos, omega_l, 'green', 'LineWidth', 2);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute and plot the analytical estimate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
omega_est = zeros(1,m);
for i = 1:m
    K_vco = K_vcos_semistable_cycle(i);
    xi = (k*tau_2*K_vco + 1)  /(2*sqrt(K_vco*(tau_1 + tau_2)));
    eta = (k*tau_2*K_vco - mu)/(2*sqrt(K_vco*(tau_1 + tau_2)));
    rho = sqrt(abs(xi^2 - k));
    kappa = sqrt(eta^2 + k*mu);
    if(K_vco < K_vco_fn)
        r(i) = (k*tau_2*K_vco + 2*xi*sqrt(tau_2*K_vco) + 1) /...
            (k*tau_2*K_vco - 2*xi*sqrt(tau_2*K_vco) + 1) *...
    exp(2*xi/rho*(atan((1 - k*tau_2*K_vco) / (2*rho*sqrt(tau_2*K_vco))) + pi/2));
    elseif(K_vco == K_vco_fn)
        r(i) = ((k*sqrt(tau_2*K_vco) + sqrt(k)) ./ (k*sqrt(tau_2*K_vco) + sqrt(k)) *...
    exp(2*sqrt(k*tau_2*K_vco) / (k*tau_2*K_vco - 1)))^2;
    else
        r(i) = ((k*sqrt(tau_2*K_vco) + xi)^2 - rho^2) / ((k*sqrt(tau_2*K_vco) - xi)^2 - rho^2) *...
    (((k*sqrt(tau_2*K_vco) + rho)^2 - xi^2) / ((k*sqrt(tau_2*K_vco) - rho)^2 - xi^2))^(xi/rho);
    end
    omega_est(i) = (sqrt(r(i))-1)/(sqrt(r(i))+1)*K_vco;
end
    
    
semilogx(K_vcos_semistable_cycle, omega_est./K_vcos_semistable_cycle, 'red', 'LineWidth', 2);


% Formatting the plot
set(gca, 'FontSize', 15);
xlabel('\textbf{$K_{\rm vco}$}', 'Interpreter', 'latex', 'fontsize', 20);
ylabel('\textbf{$\frac{\omega_e^{\rm free}}{K_{\rm vco}}$}', 'Interpreter', 'latex', 'fontsize', 20);
axis([min(K_vcos), max(K_vcos), 0, 1.1]); % Set axis limits
 


 
  
        
 
