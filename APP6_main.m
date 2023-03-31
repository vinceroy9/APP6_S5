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


%% Identification

figure(1)
plot(t,acc_mes,'kx')
xlabel('Temps [s]')
ylabel('Accélération (D_{aéro} / m) en [m/s^2]')
legend('Points mesurés','Location','NorthWest')
title('Mesures accélérométriques de la NASA avec \gamma = -90 deg')





% Intégrons les données pour obtenir v(h). On cherche donc deux set de
% données qui donne des v et des h pour des temps identiques qu'on met
% ensemble par la suite. a->v par Trapeze et v->h par Simpson. Int

% 1. Trapèze:
v_trap(1) = 0;  %Initialisation du vecteur
dt = t(2)-t(1);
%dt = (t(end) - t(1))/length(t)
for n = 2:length(acc_mes)
    v_trap(n) = dt/2 * (acc_mes(1) + acc_mes(n) + 2*sum(acc_mes(2:n-1)));
end
v_trap = -v_trap + v_ini_nasa; % Pour ajouter la vitesse initiale à t = 0;

% Erreur sur trapeze
fpa_v = (acc_mes(2)-acc_mes(1))/dt; 
fpb_v = (acc_mes(end)-acc_mes(end-1))/dt; 
err_v = abs(dt^2/12 * (fpb_v - fpa_v));

% 2. Simpson:
h_simp(1) = 0;  %Initialisation du vecteur
t_simp(1) = t(1);
for n = 3:2:length(v_trap)
    h_simp((n-1)/2+1) = dt/3 * (-v_trap(1) + -v_trap(n) + 4*sum(-v_trap(2:2:n-1)) + 2*sum(-v_trap(3:2:n-1)));
    t_simp((n-1)/2+1) = t(n);
end
h_simp = h_simp + h_ini_nasa; % Pour ajouter la hauteur initiale à t = 0;

% Erreur sur Simpson
fpppa_h = ( -v_trap(4) - -3*v_trap(3) + -3*v_trap(2) - -v_trap(1) )/dt^3; 
fpppb_h = ( -v_trap(end) - -3*v_trap(end-1) + -3*v_trap(end-2) - -v_trap(end-3) )/dt^3;
err_h = abs(dt^4/180 * (fpppb_h - fpppa_h));

% Les paires de données sont donc v_trap(1:2:end) et h_simp pour t_simp (ou
% t(1:2:end)

% 3. Approximation à 2 paramètres:
Pdyn = acc_mes*m/(S*CD0);

Y = log(2*Pdyn(1:2:end)./(v_trap(1:2:end)'.^2));
X = h_simp';

b_m = inv([length(X),sum(X);sum(X),sum(X.^2)])*[sum(Y);sum(Y.*X)];

hs = -1/b_m(2)      % hs = 11100 [m] ou 1.11e4 [km] selon google
rho0 = exp(b_m(1))  % rho0 = 0.020 kg/m^3 selon google

