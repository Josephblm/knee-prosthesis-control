%% KNEE PROSTHESIS CONTROL SIMULATION
% Joseph Belamich - EPFL Exchange Student at HKUST
% LQR Control of a Bionic Knee Prosthesis
clear all; close all; clc;

%% Parametres physiques
I = 0.5;    % Moment d'inertie (kg.m2)
b = 2.0;    % Amortissement (N.m.s/rad)
k = 10.0;   % Raideur (N.m/rad)

%% Espace d'etat
A = [0,    1;
    -k/I, -b/I];
B = [0; 1/I];

%% Controleur LQR
Q = diag([2000, 50]);
R = 0.01;
K = lqr(A, B, Q, R);

%% Trajectoire de reference - cycle de marche naturel
theta_ref_func = @(t) ...
    (t < 0.5) .* 0 + ...
    (t >= 0.5 & t < 0.75) .* (60 * (1 - cos(pi*(t-0.5)/0.25))/2) + ...
    (t >= 0.75 & t <= 1.0) .* (60 * (1 + cos(pi*(t-0.75)/0.25))/2);

theta_dot_ref_func = @(t) ...
    (t >= 0.5 & t < 0.75) .* (60*pi/(0.25*2) * sin(pi*(t-0.5)/0.25)) + ...
    (t >= 0.75 & t <= 1.0) .* (-60*pi/(0.25*2) * sin(pi*(t-0.75)/0.25));

%% Simulation avec ode45
dynamics = @(t, x) A*x + B*(-K*(x - ...
    [theta_ref_func(t)*pi/180; theta_dot_ref_func(t)*pi/180]));

x0 = [0; 0];
[t, x] = ode45(dynamics, [0 1.0], x0);

%% Trajectoire de reference sur les memes instants
theta_ref = theta_ref_func(t);

%% Visualisation
figure;

subplot(3,1,1);
plot(t, theta_ref, 'g--', 'LineWidth', 2.5); hold on;
plot(t, x(:,1)*180/pi, 'b', 'LineWidth', 2);
xlabel('Temps (s)');
ylabel('Angle (deg)');
title('Trajectoire du genou avec controleur LQR');
legend('Reference (marche naturelle)', 'Prothese controlee', 'Location', 'best');
grid on;
ylim([-5 70]);

subplot(3,1,2);
plot(t, x(:,2)*180/pi, 'r', 'LineWidth', 2);
xlabel('Temps (s)');
ylabel('Vitesse angulaire (deg/s)');
title('Vitesse angulaire');
grid on;

subplot(3,1,3);
tau = zeros(length(t),1);
for i = 1:length(t)
    e = x(i,:)' - [theta_ref_func(t(i))*pi/180; ...
        theta_dot_ref_func(t(i))*pi/180];
    tau(i) = -K*e;
end
plot(t, tau, 'm', 'LineWidth', 2);
xlabel('Temps (s)');
ylabel('Couple (N.m)');
title('Commande du moteur');
grid on;

%% COMPARAISON AVEC / SANS CONTROLEUR
% Simulation sans controleur (genou libre)
x0_free = [0; 0];
dynamics_free = @(t, x) A*x;  % pas de commande

% On part du meme point de depart que la reference a t=0.5s
x0_free = [theta_ref_func(0.5)*pi/180; 0];
[t_free, x_free] = ode45(dynamics_free, [0.5 1.0], x0_free);

% Trajectoire reference sur ces instants
theta_ref_free = theta_ref_func(t_free);

%% Nouveau graphique comparaison
figure;
plot(t, theta_ref, 'g--', 'LineWidth', 2.5); hold on;
plot(t, x(:,1)*180/pi, 'b', 'LineWidth', 2.5);
plot(t_free, x_free(:,1)*180/pi, 'r', 'LineWidth', 2);
xlabel('Temps (s)');
ylabel('Angle du genou (deg)');
title('Comparaison : avec et sans controleur LQR');
legend('Reference (marche naturelle)', ...
       'Prothese avec LQR', ...
       'Prothese sans controleur', ...
       'Location', 'best');
grid on;
ylim([-10 70]);


%% OPTION B - PLUSIEURS CYCLES DE MARCHE
% On simule 4 cycles consecutifs (4 pas)
n_cycles = 4;
T_cycle = 1.0;  % duree d'un cycle

% Trajectoire de reference sur plusieurs cycles
theta_ref_multi = @(t) theta_ref_func(mod(t, T_cycle));
theta_dot_ref_multi = @(t) theta_dot_ref_func(mod(t, T_cycle));

% Simulation sur 4 cycles
dynamics_multi = @(t, x) A*x + B*(-K*(x - ...
    [theta_ref_multi(t)*pi/180; ...
     theta_dot_ref_multi(t)*pi/180]));

[t_multi, x_multi] = ode45(dynamics_multi, [0 n_cycles*T_cycle], [0;0]);

% Reference sur ces instants
theta_ref_m = arrayfun(theta_ref_multi, t_multi);

% Visualisation
figure;
plot(t_multi, theta_ref_m, 'g--', 'LineWidth', 2); hold on;
plot(t_multi, x_multi(:,1)*180/pi, 'b', 'LineWidth', 2);
xlabel('Temps (s)');
ylabel('Angle du genou (deg)');
title('4 cycles de marche consecutifs avec controleur LQR');
legend('Reference', 'Prothese controlee', 'Location', 'best');
grid on;
ylim([-5 70]);

% Ajouter lignes verticales pour separer les cycles
for i = 1:n_cycles-1
    xline(i*T_cycle, 'k--', 'LineWidth', 1, ...
          'Label', sprintf('Cycle %d', i));
end

%% OPTION C - ROBUSTESSE AU BRUIT
% Dans la vraie vie les capteurs ont du bruit
% On simule un bruit gaussien sur la mesure d'angle

noise_level = 0.5;  % bruit en degres (ecart type)

% Simulation avec bruit sur les capteurs
dt = 0.001;
tspan_n = 0:dt:1.0;
N = length(tspan_n);

x_noisy = zeros(N, 2);
x_noisy(1,:) = [0; 0];
tau_noisy = zeros(N, 1);

for i = 1:N-1
    % Mesure bruitee (ce que le capteur lit reellement)
    noise = noise_level * pi/180 * randn(2,1);
    x_measured = x_noisy(i,:)' + noise;
    
    % Reference a cet instant
    x_des = [theta_ref_func(tspan_n(i))*pi/180; ...
             theta_dot_ref_func(tspan_n(i))*pi/180];
    
    % Commande basee sur la mesure bruitee
    e = x_measured - x_des;
    tau = -K * e;
    tau_noisy(i) = tau;
    
    % Evolution reelle du systeme
    x_dot = A * x_noisy(i,:)' + B * tau;
    x_noisy(i+1,:) = x_noisy(i,:)' + x_dot * dt;
end

% Reference
theta_ref_n = theta_ref_func(tspan_n);

% Visualisation
figure;
plot(tspan_n, theta_ref_n, 'g--', 'LineWidth', 2.5); hold on;
plot(t, x(:,1)*180/pi, 'b', 'LineWidth', 2);
plot(tspan_n, x_noisy(:,1)*180/pi, 'r', 'LineWidth', 1.5);
xlabel('Temps (s)');
ylabel('Angle du genou (deg)');
title('Robustesse du controleur LQR face au bruit des capteurs');
legend('Reference', 'LQR sans bruit', ...
       'LQR avec bruit capteurs', 'Location', 'best');
grid on;
ylim([-5 70]);

%% OPTION D - COMPARAISON LQR vs PID
% Le PID est le controleur classique utilise en industrie
% On va montrer que le LQR est plus performant

% Parametres PID (tuning manuel)
Kp = 80;    % gain proportionnel
Ki = 5;     % gain integral
Kd = 3;     % gain derivateur

% Simulation PID
dt = 0.001;
tspan_pid = 0:dt:1.0;
N = length(tspan_pid);

x_pid = zeros(N, 2);
x_pid(1,:) = [0; 0];
integral_error = 0;

for i = 1:N-1
    % Reference
    theta_des = theta_ref_func(tspan_pid(i)) * pi/180;
    
    % Erreur
    error = x_pid(i,1) - theta_des;
    integral_error = integral_error + error * dt;
    
    if i > 1
        deriv_error = (x_pid(i,1) - x_pid(i-1,1)) / dt;
    else
        deriv_error = 0;
    end
    
    % Commande PID
    tau_pid = -(Kp*error + Ki*integral_error + Kd*deriv_error);
    
    % Evolution du systeme
    x_dot = A * x_pid(i,:)' + B * tau_pid;
    x_pid(i+1,:) = x_pid(i,:)' + x_dot * dt;
end

% Reference
theta_ref_pid = theta_ref_func(tspan_pid);

% Visualisation
figure;
plot(tspan_pid, theta_ref_pid, 'g--', 'LineWidth', 2.5); hold on;
plot(t, x(:,1)*180/pi, 'b', 'LineWidth', 2.5);
plot(tspan_pid, x_pid(:,1)*180/pi, 'r', 'LineWidth', 2);
xlabel('Temps (s)');
ylabel('Angle du genou (deg)');
title('Comparaison LQR vs PID');
legend('Reference', 'LQR (optimal)', 'PID (classique)', 'Location', 'best');
grid on;
ylim([-10 70]);

%% OPTION E - RAIDEUR VARIABLE
% Hugh Herr 2024 : prothese a raideur variable selon la phase de marche
% Phase appui : raideur elevee (genou rigide pour supporter le poids)
% Phase oscillation : raideur faible (genou souple pour swing naturel)

% Fonction de raideur variable
k_func = @(t) ...
    (mod(t,1.0) < 0.5) .* 20.0 + ...   % appui : k=20 (rigide)
    (mod(t,1.0) >= 0.5) .* 3.0;         % oscillation : k=3 (souple)

% Simulation avec raideur variable sur 2 cycles
dt = 0.001;
tspan_vs = 0:dt:2.0;
N = length(tspan_vs);

x_vs = zeros(N, 2);
x_vs(1,:) = [0; 0];
tau_vs = zeros(N, 1);
k_history = zeros(N, 1);

for i = 1:N-1
    t_i = tspan_vs(i);
    
    % Raideur a cet instant
    k_i = k_func(t_i);
    k_history(i) = k_i;
    
    % Matrices avec raideur variable
    A_i = [0, 1; -k_i/I, -b/I];
    K_i = lqr(A_i, B, Q, R);
    
    % Reference
    x_des = [theta_ref_multi(t_i)*pi/180; ...
             theta_dot_ref_multi(t_i)*pi/180];
    
    % Commande
    e = x_vs(i,:)' - x_des;
    tau_i = -K_i * e;
    tau_vs(i) = tau_i;
    
    % Evolution
    x_dot = A_i * x_vs(i,:)' + B * tau_i;
    x_vs(i+1,:) = x_vs(i,:)' + x_dot * dt;
end

% Reference sur 2 cycles
theta_ref_vs = arrayfun(theta_ref_multi, tspan_vs);

% Visualisation
figure;

subplot(3,1,1);
plot(tspan_vs, theta_ref_vs, 'g--', 'LineWidth', 2); hold on;
plot(tspan_vs, x_vs(:,1)*180/pi, 'b', 'LineWidth', 2);
ylabel('Angle (deg)');
title('Prothese a raideur variable (Hugh Herr 2024)');
legend('Reference', 'Prothese', 'Location', 'best');
grid on; ylim([-5 70]);

subplot(3,1,2);
plot(tspan_vs, k_history, 'm', 'LineWidth', 2);
ylabel('Raideur k (N.m/rad)');
title('Raideur variable selon la phase de marche');
grid on;

subplot(3,1,3);
plot(tspan_vs, tau_vs, 'r', 'LineWidth', 2);
xlabel('Temps (s)');
ylabel('Couple (N.m)');
title('Commande du moteur avec raideur variable');
grid on;

