clc
clear all
close all

load Accelero_Data_from_NASA

% ======== Paramètres ========  %
% NOTE : Mettre tout les angles en RADIANS dans les équations

m = 50;                 % kg            Masse de la capsule)
J = 1.5;                % kg-m^2        Inertie de la capsule)
R_mars = 3397e3;        % m             Rayon de mars
mu_mars = 42830e9;      % m^3/s^2       Paramètre grav. de Mars
%rho0 = ???             % ???           Densité atmosphérique à h0
%hs = ???               % ???           Facteur d'échelle de la densité
S = 0.8;                % m^2           Surface aéro. de la capsule
d = 0.05;               % m             Dimension aéro. de la capsule
CD0 = 1.2;              %               Coefficient de la trainée
CLa = 0.8;              %               Coefficient de la portance
CMa = -0.07;            %               Coefficient de couple
CMq = -0.05;            %               Coefficient d'amortissement
CMd = 0.1;              %               Coefficient de volet aéro.

% ======== Condition intiales ======== %

v_ini = 6100;           % m/s
gamma_ini = -20.5;      % deg       !!!
h_ini = 120000;         % m
s_ini = 0.0;            % deg       !!!
theta_ini = -80;        % deg       !!!
q_ini = 0.0;            % deg/s     !!!

% ======== Condition finales désirées ======== %

v_fin = 250;            % m/s       OU
%v_fin = 300;           % m/s
h_fin = 10000;          % m

% ======== Contraintes ======== %

% Delta_t_lim (D_aero > 2650 N) <= 40 s
% P_dyn_max <= 9500 N/m^2
% │theta_cmd│ < 60 deg


%% Résultat de la NASA

% ======== Condition intiales ======== %

v_ini_nasa = 6100;          % m/s
gamma_ini_nasa = -90;       % deg       
h_ini_nasa = 120000;        % m
s_ini_nasa = 0.0;           % deg       
theta_ini_nasa = -90;       % deg       
q_ini_nasa = 0.0;           % deg/s     


% ======== Bruit de mesure et Paramètre Balistique ======== %
sigma_n = 0.035; % ??? car sigma_n((a_mes-a_approx)/a_mes) donc mult ???



figure(1)
plot(t,acc_mes,'kx')
xlabel('Temps [s]')
ylabel('Accélération (D_{aéro} / m) en [m/s^2]')
legend('Points mesurés','Location','NorthWest')
title('Mesures accélérométriques de la NASA avec \gamma = -90 deg')







