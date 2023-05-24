close all;
clc;
clear all;

k = я/pi;
mu = pi*k - 1;


% array_a_sawtooth = [0, 0.001, 0.01, 0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 0.9];
% array_a = [0, 0.005, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.7, 0.9];
a_array = [0.9, 0.7, 0.5, 0.4, 0.3, 0.2, 0.1, 0.05, 0.01, 0.001, 0];
% a_array = [0.2];



%btw, standard engineering parameters tau_1 = 0.0448, tau_2 = 0.0185 correspond to a = 0.2923

x_array = [0.0001:0.01:10 11:0.1:100 110:1:1000 1100:10:10000];
% x_array = [506];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Hold-in %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h = figure();
semilogx([min(x_array) max(x_array)], [1 1], 'm', 'LineWidth', 3);
hold on;
grid on;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Pull-in %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:length(a_array)
    a = a_array(i);

xi = (a*k*x_array + 1)./(2*sqrt(x_array));
eta = (a*k*x_array - mu)./(2*sqrt(x_array));
rho = sqrt(abs(xi.^2 - k));
kappa = sqrt(eta.^2 + mu*k);



if a==0
    x_sep = x_array(x_array(:) > min(x_array));
    y_sep = omega_sep_for_diagram(a, x_sep, k);
%     semilogx(x_minus, 1, 'x','LineWidth', 2, 'Color', 'black');
    semilogx(x_sep, y_sep, 'black', 'LineWidth', 1);
    y_l = omega_l_for_diagram(a,x_array,k);
    semilogx(x_array(x_array(:) >= min(x_array)), y_l(x_array(:) >= min(x_array)), 'blue', 'LineWidth', 1);
else
    x_ht = 1/(k*(2 - a + 2*sqrt(1 - a)));
    x_pt = mu/(a*k);
    x_nf = 1/(k*(2 - a - 2*sqrt(1 - a)));

    %         y_lyapunov = (1 + sqrt(a) - sqrt((1 + sqrt(a))^2 - 4*a))/(2*sqrt(a));
%     semilogx([min(x_array) max(x_array)], [y_lyapunov y_lyapunov], 'red', 'LineWidth', 1);
    
    semilogx([min(x_array) x_ht], [1 1], 'black', 'LineWidth', 1);
    
    x_sep = x_array(x_array(:) >= x_ht & x_array(:) <= x_pt);
    y_sep = omega_sep_for_diagram(a, x_sep, k);
    semilogx(x_sep, y_sep, 'black', 'LineWidth', 1);
    semilogx(x_pt, omega_sep_for_diagram(a, x_pt, k), 'x','LineWidth', 2, 'Color', 'black');
    
    semilogx(x_nf, omega_ss_for_diagram(a, x_nf, k), 'x','LineWidth', 2, 'Color', 'red');


x_sep = x_array(x_array(:) >= x_pt);
    y_sep = omega_sep_for_diagram(a, x_sep, k);
    semilogx(x_sep, y_sep, 'black--', 'LineWidth', 1);
    semilogx(x_nf, omega_sep_for_diagram(a, x_nf, k), 'x','LineWidth', 2, 'Color', 'blue');


    x_ss = x_array(x_array(:) > max(x_pt, x_ht));
    y_ss = omega_ss_for_diagram(a, x_ss, k);
    semilogx(x_ss, y_ss, 'black', 'LineWidth', 1);


   





    %%%Safonov's asymptotic formula
    fcn = @(theta) (theta^2/(sinh(theta))^2 - (1-a));
    start=eps;
    finish = 1/eps;
    theta_0 = fzero(fcn, [start, finish]);
    y_as = (sinh(theta_0)*cosh(theta_0)-theta_0)/(sinh(theta_0))^2
    semilogx([1 max(x_array)], [y_as y_as], 'red--', 'LineWidth', 1);

   


    
    %%%analytic estimate
    r_1 = ((k*sqrt(a*x_array) + eta).^2 - kappa.^2) ./ ((k*sqrt(a*x_array) - eta).^2 - kappa.^2) .*...
    (((k*sqrt(a*x_array) + kappa).^2 - eta.^2) ./ ((k*sqrt(a*x_array) - kappa).^2 - eta.^2)).^(eta./kappa);

    x_21 = x_array(x_array(:) < x_nf);
    xi_21 = xi(x_array(:) < x_nf);
    rho_21 = rho(x_array(:) < x_nf);
    r_21 = (k*a*x_21 + 2*xi_21.*sqrt(a*x_21) + 1) ./ (k*a*x_21 - 2*xi_21.*sqrt(a*x_21) + 1) .*...
    exp(2*xi_21./rho_21.*(atan((1 - k*a*x_21) ./ (2*rho_21.*sqrt(a*x_21))) + pi/2));

    x_22 = x_array(x_array(:) > x_nf);
    xi_22 = xi(x_array(:) > x_nf);
    rho_22 = rho(x_array(:) > x_nf);
    r_22 = ((k*sqrt(a*x_22) + xi_22).^2 - rho_22.^2) ./...
        ((k*sqrt(a*x_22) - xi_22).^2 - rho_22.^2) .*...
    (((k*sqrt(a*x_22) + rho_22).^2 - xi_22.^2) ./ ...
    ((k*sqrt(a*x_22) - rho_22).^2 - xi_22.^2)).^(xi_22./rho_22);
    
    r_2 = [r_21, r_22];
    r = max(r_1, r_2);
    y = (sqrt(r) - 1)./(sqrt(r) + 1);
%         semilogx(x_array, y, 'g', 'LineWidth', 1);
         


%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Lock-in %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     y_l = omega_l_for_diagram(a,x_array,k);
%     semilogx(x_array(x_array(:) >= min(x_array)), y_l(x_array(:) >= min(x_array)), 'blue', 'LineWidth', 1);
end    
end


set(gca,'FontSize', 15);

xlabel('\textbf{$(\tau_1 + \tau_2)K_{\rm vco}$}','Interpreter','latex', 'fontsize', 20);  
ylabel('\textbf{$\frac{\omega_l}{K_{\rm vco}}$}','Interpreter','latex', 'fontsize', 20); 
% legend('\omega_h', '\omega_{\rm sep}','\omega_{\rm semi}', '\omega_l', 'Location','best');
axis([0.09 max(x_array) 0 1.01]);
xticks([0.01 0.1 1 10 100 1000 5000]);
xticklabels({'10^{-2}','10^{-1}','10^0','10^1','10^2','10^3','5\cdot10^3'});



