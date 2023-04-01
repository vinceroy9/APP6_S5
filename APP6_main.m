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

v_fin1 = 250;           % m/s       
v_fin2 = 300;           % m/s
h_fin = 10000;          % m

% ======== Contraintes ======== %

Delta_t_lim = 40;        % s
P_dyn_max = 9500;        % N/m^2
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
sigma_n = 0.035;            % ??? car sigma_n((a_mes-a_approx)/a_mes) donc mult ???
B = CD0*S/m;                % m^2/kg

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

figure
plot(t, v_trap)

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

% 4.  Calcul de l'erreur RMS et R2
% RMS
a_approx = (S*CD0)/m * 1/2*rho0*exp(-h_simp/hs).*v_trap(1:2:end).^2;
E = sum((a_approx-acc_mes(1:2:end)').^2);
RMS_a = sqrt(1/length(acc_mes(1:2:end)) * E)   % NOTE : RMS de l'acceleration

% NOTE : R2 TOUJOURS calculée dans le domaine linéaire (Y,X)
Y_calc = log(rho0*exp(-h_simp/hs));
R2 = sum((Y_calc-mean(Y)).^2)/sum((Y-mean(Y)).^2);


% 5. Comparaison avec bruit de l'accelerometre + RMS absolue (a_mes-a_approx)
% et (a_mes-a_approx)/a_mes

% ????????????


%% Validation de la RAA
h_raa = linspace(h_simp(1),h_simp(end),1e4);

rho_raa = 1/2 * rho0 * exp(-h_raa/hs);
rho_raa_ini = 1/2 * rho0 * exp(-h_ini_nasa/hs);
v_raa = v_ini_nasa * exp(1/2 * B * hs * (rho_raa - rho_raa_ini)/sind(gamma_ini_nasa));

figure(2)
plot(h_simp,v_trap(1:2:end),'xb')
hold on
plot(h_raa,v_raa,'r')
hold off
xlabel('Altidude h [m]')
ylabel('Vitesse v [m/s]')
legend('v(h)','RAA','Location','NorthWest')

%% Vérification des limites structurelle de la capsule

% Exprimons la courbe Pdyn avec RAA+gravité et appliquons Newton-Rhapson
% sur celui-ci. Note: Ce Pdyn sera en fonction de h et le delta_t_lim sera
% fait a partir du delta_h trouver par newton-rhapson.


% 0. Calcu de delta_v_aero
% delta_v_aero1 = v_fin1 - sqrt(v_raa.^2 + 2*mu_mars*(1/(h_fin+R_mars) - 1./(h_raa+R_mars)));
% delta_v_aero2 = v_fin2 - sqrt(v_raa.^2 + 2*mu_mars*(1/(h_fin+R_mars) - 1./(h_raa+R_mars)));

delta_v_aero1 = v_fin1 - sqrt(v_ini^2 + 2*mu_mars*(1/(h_fin+R_mars) - 1/(h_ini+R_mars)));
delta_v_aero2 = v_fin2 - sqrt(v_ini^2 + 2*mu_mars*(1/(h_fin+R_mars) - 1/(h_ini+R_mars)));

% 1. Calcul de l'angle gamma_ref pour les deux vitesse terminales

% Rho final utilisé dans les calculs de gamma ref
rho_fin = rho0*exp(-h_fin/hs);

gamma_ref1 = asind(0.5*B*hs*(rho_fin - rho_raa_ini)/(log(1 + delta_v_aero1/v_ini))); 
gamma_ref2 = asind(0.5*B*hs*(rho_fin - rho_raa_ini)/(log(1 + delta_v_aero2/v_ini))); 

% Ça c'est le vecteur avec tous les gamma ref calculés mais je sais pas
% lequel je dois prendre tho (PER PHILIPPE, PRENDRE LE INITIAL)
% gamma_ref1 = asind(0.5*B*hs*(rho_fin - rho_raa)./(log(1 + delta_v_aero1./v_raa)));
% gamma_ref2 = asind(0.5*B*hs*(rho_fin - rho_raa)./(log(1 + delta_v_aero2./v_raa)));

% Je refais le vecteur de v avec la raa
v_raa_ref1 =  v_ini_nasa * exp(1/2 * B * hs * (rho_raa - rho_raa_ini)./sind(gamma_ref1));
v_raa_ref2 =  v_ini_nasa * exp(1/2 * B * hs * (rho_raa - rho_raa_ini)./sind(gamma_ref2));

Pdyn_ref1 = 1/2* rho_raa .* v_raa_ref1.^2;
Pdyn_max1 = max(Pdyn_ref1);

Pdyn_ref2 = 1/2* rho_raa .* v_raa_ref2.^2;
Pdyn_max2 = max(Pdyn_ref2); % ON EST GOOD C'EST EN BAS DE 9500

% 2. Newton-Raphson 

% Je suis rendu ici :') (2023-04-01 18:03)






