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
D_aero_max = 2650;       % N (pas dépasser ça plus que 40s)
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

% figure
% plot(t, v_trap)

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

rho_raa = rho0 * exp(-h_raa/hs);
rho_raa_ini = rho0 * exp(-h_ini_nasa/hs);
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
delta_v_aero1 = v_fin1 - sqrt(v_ini^2 + 2*mu_mars*(1/(h_fin+R_mars) - 1/(h_ini+R_mars)));
delta_v_aero2 = v_fin2 - sqrt(v_ini^2 + 2*mu_mars*(1/(h_fin+R_mars) - 1/(h_ini+R_mars)));

% 1. Calcul de l'angle gamma_ref pour les deux vitesse terminales
% Rho final utilisé dans les calculs de gamma ref
rho_fin = rho0*exp(-h_fin/hs);

% Calculs gamma refs
gamma_ref1 = asind(0.5*B*hs*(rho_fin - rho_raa_ini)/(log(1 + delta_v_aero1/v_ini))); 
gamma_ref2 = asind(0.5*B*hs*(rho_fin - rho_raa_ini)/(log(1 + delta_v_aero2/v_ini))); 

% Je refais les vecteurs de v avec la raa
v_raa_ref1 =  v_ini_nasa * exp((1/2) * B * hs * (rho_raa - rho_raa_ini)./sind(gamma_ref1));
v_raa_ref2 =  v_ini_nasa * exp((1/2) * B * hs * (rho_raa - rho_raa_ini)./sind(gamma_ref2));

% Pour vfin1 = 250 m/s
P_dyn_ref1 = (1/2)* rho_raa .* v_raa_ref1.^2;
D_aero_ref1 = P_dyn_ref1*S*CD0;
P_dyn_max1 = max(P_dyn_ref1); % ON EST GOOD C'EST EN BAS DE 9500
% F_ref1 = D_aero_ref1 - D_aero_max;

% Pour vfin2 = 300 m/s
P_dyn_ref2 = (1/2)* rho_raa .* v_raa_ref2.^2;
D_aero_ref2 = P_dyn_ref2*S*CD0;
P_dyn_max2 = max(P_dyn_ref2); % ON EST GOOD C'EST EN BAS DE 9500

% Affichage des Daero en fonction de l'altitde
figure
plot(h_raa, D_aero_ref1)
hold on
plot(h_raa, ones(length(h_raa), 1)*D_aero_max)
title('Daero selon altitude :  vfin = 250 m/s')
xlabel('Altitude h [m]')
ylabel('Trainee Daero [N]') 
hold off

figure
plot(h_raa, D_aero_ref2)
hold on
plot(h_raa, ones(length(h_raa), 1)*D_aero_max)
title('Daero selon altitude :  vfin = 300 m/s')
xlabel('Altitude h [m]')
ylabel('Trainee Daero [N]') 
hold off

% 2. Newton-Raphson
tol = 1e-10; % Donnée dans guide étudiant

% =======================     VFIN1 = 250 m/s   ========================= %
% ===================     c'est avec gamma_ref_1   ====================== %
% points de départ des itérations, voir graph pour comprendre
h_vfin1_1 = 19000; 
h_vfin1_2 = 40000; 

% rho est fonction de h
rho_vfin1_1 = rho0 * exp(-h_vfin1_1/hs);
rho_vfin1_2 = rho0 * exp(-h_vfin1_2/hs);

% v aussi fonction de h (Pour RAA)
v_vfin1_1 = v_ini_nasa * exp(1/2*B*hs*(rho_vfin1_1 - rho_raa_ini)/sind(gamma_ref1));
v_vfin1_2 = v_ini_nasa * exp(1/2*B*hs*(rho_vfin1_2 - rho_raa_ini)/sind(gamma_ref1));

% donc Pdyn est fonction de h
Pdyn_vfin1_1 = (1/2)*rho_vfin1_1*v_vfin1_1^2;
Pdyn_vfin1_2 = (1/2)*rho_vfin1_2*v_vfin1_2^2;

% donc Daero aussi
D_aero_vfin1_1 = Pdyn_vfin1_1*S*CD0;
D_aero_vfin1_2 = Pdyn_vfin1_2*S*CD0;

% On veut les valeurs qui croisent Daero max
F_vfin1_1 = D_aero_vfin1_1 - D_aero_max;
F_vfin1_2 = D_aero_vfin1_2 - D_aero_max;

% Dérivée de Daero par rapport à h (un peu dégueulasse)
% Je l'ai split en parties
diff_rho_vfin1_1 = (-rho0/hs)*exp(-h_vfin1_1/hs);
diff_v2_vfin1_1 = ((-v_ini_nasa^2*B*hs*rho0)/(sind(gamma_ref1)*hs)) * exp(-h_vfin1_1/hs) * exp(B*hs*(rho_vfin1_1-rho_raa_ini)/sind(gamma_ref1));

diff_rho_vfin1_2 = (-rho0/hs)*exp(-h_vfin1_2/hs);
diff_v2_vfin1_2 = ((-v_ini_nasa^2*B*hs*rho0)/(sind(gamma_ref1)*hs)) * exp(-h_vfin1_2/hs) * exp(B*hs*(rho_vfin1_2-rho_raa_ini)/sind(gamma_ref1));

% Voici le résultat
D_vfin1_1 = (S*CD0/2)*(diff_rho_vfin1_1*v_vfin1_1^2 + diff_v2_vfin1_1*rho_vfin1_1);
D_vfin1_2 = (S*CD0/2)*(diff_rho_vfin1_1*v_vfin1_2^2 + diff_v2_vfin1_2*rho_vfin1_2);

% Itérations Newton Raphson
% Premier h 
count1 = 0;
while abs(F_vfin1_1) > tol
    % Nouveau h selon Newton-Raphson
    h_vfin1_1 = h_vfin1_1 - F_vfin1_1/D_vfin1_1;
    
    % Recalcul avec ce h -> rho -> v -> Pdyn -> Daero -> F ->dérivée Daero
    rho_vfin1_1 = rho0 * exp(-h_vfin1_1/hs);
    v_vfin1_1 = v_ini_nasa * exp(1/2*B*hs*(rho_vfin1_1 - rho_raa_ini)/sind(gamma_ref1));
    Pdyn_vfin1_1 = 1/2*rho_vfin1_1*v_vfin1_1^2;
    D_aero_vfin1_1 = Pdyn_vfin1_1*S*CD0;
    F_vfin1_1 = D_aero_vfin1_1 - D_aero_max;
    
    % Partie utilisée dans dérivée
    diff_rho_vfin1_1 = (-rho0/hs)*exp(-h_vfin1_1/hs);
    diff_v2_vfin1_1 = ((-v_ini_nasa^2*B*hs*rho0)/(sind(gamma_ref1)*hs)) * exp(-h_vfin1_1/hs) * exp(B*hs*(rho_vfin1_1-rho_raa_ini)/sind(gamma_ref1));
    
    D_vfin1_1 = (S*CD0/2)*(diff_rho_vfin1_1*v_vfin1_1^2 + diff_v2_vfin1_1*rho_vfin1_1);
    
    % Compteur pour éviter boucle infinie
    count1 = count1 + 1;
    if count1 > 500
        break;
    end
end

% Deuxieme h 
count2 = 0;
while abs(F_vfin1_2) > tol
    % Nouveau h selon Newton-Raphson
    h_vfin1_2 = h_vfin1_2 - F_vfin1_2/D_vfin1_2;
    
    % Recalcul avec ce h -> rho -> v -> Pdyn -> Daero -> F ->dérivée Daero
    rho_vfin1_2 = rho0 * exp(-h_vfin1_2/hs);
    v_vfin1_2 = v_ini_nasa * exp(1/2*B*hs*(rho_vfin1_2 - rho_raa_ini)/sind(gamma_ref1));
    Pdyn_vfin1_2 = 1/2*rho_vfin1_2*v_vfin1_2^2;
    D_aero_vfin1_2 = Pdyn_vfin1_2*S*CD0;
    F_vfin1_2 = D_aero_vfin1_2 - D_aero_max;
    
    % Partie utilisée dans dérivée
    diff_rho_vfin1_2 = (-rho0/hs)*exp(-h_vfin1_2/hs);
    diff_v2_vfin1_2 = ((-v_ini_nasa^2*B*hs*rho0)/(sind(gamma_ref1)*hs)) * exp(-h_vfin1_2/hs) * exp(B*hs*(rho_vfin1_2-rho_raa_ini)/sind(gamma_ref1));
    
    D_vfin1_2 = (S*CD0/2)*(diff_rho_vfin1_2*v_vfin1_2^2 + diff_v2_vfin1_2*rho_vfin1_2);
    
    % Compteur pour éviter boucle infinie
    count2 = count2 + 1;
    if count2 > 500
        break;
    end
end

% Vérification de la dérivée
% for i = 1 : length(F_ref1)-1
%     diff_approx(i) = (F_ref1(i+1) - F_ref1(i))/(h_raa(2)-h_raa(1));
%     
%     diff_rho = (-rho0/hs)*exp(-h_raa(i)/hs);
%     diff_v2 = ((-B*hs*rho0*v_ini_nasa^2)/(sind(gamma_ref1)*hs)) * exp(-h_raa(i)/hs) * exp(B*hs*(rho_raa(i)-rho_raa_ini)/sind(gamma_ref1));
%     
%     D_NR1(i) = (S*CD0/2)*(diff_rho*v_raa_ref1(i)^2 + diff_v2*rho_raa(i));
% end

% Graph dérivée et approximation 
% figure
% plot(h_raa(1:end-1), D_NR1)
% hold on
% plot(h_raa(1:end-1), diff_approx)
% legend('Dérivée trouvée analytiquement', 'Dérivée approx.')
% title('Comparaison des dérivées')

% =======================     VFIN1 = 300 m/s   ========================= %
% ===================     c'est avec gamma_ref_2   ====================== %
% points de départ des itérations, voir graph pour comprendre
h_vfin2_1 = 19000; 
h_vfin2_2 = 40000; 

% rho est fonction de h
rho_vfin2_1 = rho0 * exp(-h_vfin2_1/hs);
rho_vfin2_2 = rho0 * exp(-h_vfin2_2/hs);

% v aussi fonction de h (Pour RAA)
v_vfin2_1 = v_ini_nasa * exp(1/2*B*hs*(rho_vfin2_1 - rho_raa_ini)/sind(gamma_ref2));
v_vfin2_2 = v_ini_nasa * exp(1/2*B*hs*(rho_vfin2_2 - rho_raa_ini)/sind(gamma_ref2));

% donc Pdyn est fonction de h
Pdyn_vfin2_1 = (1/2)*rho_vfin2_1*v_vfin2_1^2;
Pdyn_vfin2_2 = (1/2)*rho_vfin2_2*v_vfin2_2^2;

% donc Daero aussi
D_aero_vfin2_1 = Pdyn_vfin2_1*S*CD0;
D_aero_vfin2_2 = Pdyn_vfin2_2*S*CD0;

% On veut les valeurs qui croisent Daero max
F_vfin2_1 = D_aero_vfin2_1 - D_aero_max;
F_vfin2_2 = D_aero_vfin2_2 - D_aero_max;

% Dérivée de Daero par rapport à h (un peu dégueulasse)
% Je l'ai split en parties
diff_rho_vfin2_1 = (-rho0/hs)*exp(-h_vfin2_1/hs);
diff_v2_vfin2_1 = ((-v_ini_nasa^2*B*hs*rho0)/(sind(gamma_ref2)*hs)) * exp(-h_vfin2_1/hs) * exp(B*hs*(rho_vfin2_1-rho_raa_ini)/sind(gamma_ref2));

diff_rho_vfin2_2 = (-rho0/hs)*exp(-h_vfin2_2/hs);
diff_v2_vfin2_2 = ((-v_ini_nasa^2*B*hs*rho0)/(sind(gamma_ref2)*hs)) * exp(-h_vfin2_2/hs) * exp(B*hs*(rho_vfin2_2-rho_raa_ini)/sind(gamma_ref2));

% Voici le résultat
D_vfin2_1 = (S*CD0/2)*(diff_rho_vfin2_1*v_vfin2_1^2 + diff_v2_vfin2_1*rho_vfin2_1);
D_vfin2_2 = (S*CD0/2)*(diff_rho_vfin2_1*v_vfin2_2^2 + diff_v2_vfin2_2*rho_vfin2_2);

% Itérations Newton Raphson
% Premier h 
count3 = 0;
while abs(F_vfin2_1) > tol
    % Nouveau h selon Newton-Raphson
    h_vfin2_1 = h_vfin2_1 - F_vfin2_1/D_vfin2_1;
    
    % Recalcul avec ce h -> rho -> v -> Pdyn -> Daero -> F ->dérivée Daero
    rho_vfin2_1 = rho0 * exp(-h_vfin2_1/hs);
    v_vfin2_1 = v_ini_nasa * exp(1/2*B*hs*(rho_vfin2_1 - rho_raa_ini)/sind(gamma_ref2));
    Pdyn_vfin2_1 = 1/2*rho_vfin2_1*v_vfin2_1^2;
    D_aero_vfin2_1 = Pdyn_vfin2_1*S*CD0;
    F_vfin2_1 = D_aero_vfin2_1 - D_aero_max;
    
    % Partie utilisée dans dérivée
    diff_rho_vfin2_1 = (-rho0/hs)*exp(-h_vfin2_1/hs);
    diff_v2_vfin2_1 = ((-v_ini_nasa^2*B*hs*rho0)/(sind(gamma_ref2)*hs)) * exp(-h_vfin2_1/hs) * exp(B*hs*(rho_vfin2_1-rho_raa_ini)/sind(gamma_ref2));
    
    D_vfin2_1 = (S*CD0/2)*(diff_rho_vfin2_1*v_vfin2_1^2 + diff_v2_vfin2_1*rho_vfin2_1);
    
    % Compteur pour éviter boucle infinie
    count3 = count3 + 1;
    if count3 > 500
        break;
    end
end

% Deuxieme h 
count4 = 0;
while abs(F_vfin2_2) > tol
    % Nouveau h selon Newton-Raphson
    h_vfin2_2 = h_vfin2_2 - F_vfin2_2/D_vfin2_2;
    
    % Recalcul avec ce h -> rho -> v -> Pdyn -> Daero -> F ->dérivée Daero
    rho_vfin2_2 = rho0 * exp(-h_vfin2_2/hs);
    v_vfin2_2 = v_ini_nasa * exp(1/2*B*hs*(rho_vfin2_2 - rho_raa_ini)/sind(gamma_ref2));
    Pdyn_vfin2_2 = 1/2*rho_vfin2_2*v_vfin2_2^2;
    D_aero_vfin2_2 = Pdyn_vfin2_2*S*CD0;
    F_vfin2_2 = D_aero_vfin2_2 - D_aero_max;
    
    % Partie utilisée dans dérivée
    diff_rho_vfin2_2 = (-rho0/hs)*exp(-h_vfin2_2/hs);
    diff_v2_vfin2_2 = ((-v_ini_nasa^2*B*hs*rho0)/(sind(gamma_ref2)*hs)) * exp(-h_vfin2_2/hs) * exp(B*hs*(rho_vfin2_2-rho_raa_ini)/sind(gamma_ref2));
    
    D_vfin2_2 = (S*CD0/2)*(diff_rho_vfin2_2*v_vfin2_2^2 + diff_v2_vfin2_2*rho_vfin2_2);
    
    % Compteur pour éviter boucle infinie
    count4 = count4 + 1;
    if count4 > 500
        break;
    end
end

% NOTE: h obtenus comparés avec les graphs et ils concordent

% 3. Durée t lim 

% Création du vecteur de vraa qui comporte les éléments > v à h trouvé
v_moyen_1_vec = v_raa_ref1(v_vfin1_1 < v_raa_ref1);
v_moyen_2_vec = v_raa_ref2(v_vfin2_1 < v_raa_ref1);

% Je savais pas comment le faire en une seule ligne so je refais m chose
v_moyen_1_vec = v_moyen_1_vec(v_moyen_1_vec < v_vfin1_2);
v_moyen_2_vec = v_moyen_2_vec(v_moyen_2_vec < v_vfin2_2);

% V moyen entre les deux h
v_moyen_1 = sum(v_moyen_1_vec)/length(v_moyen_1_vec);
v_moyen_2 = sum(v_moyen_2_vec)/length(v_moyen_2_vec);

% Calculs des tlim 
delta_tlim1 = (h_vfin1_1 - h_vfin1_2)/(v_moyen_1*sind(gamma_ref1));
delta_tlim2 = (h_vfin2_1 - h_vfin2_2)/(v_moyen_2*sind(gamma_ref2));

% RÉSULTATS OBTENUS < 40s ON EST GOOD

%% Commande de la dynamique de translation
tau = 0.25; %[s]
Kp_tra = 1/tau;


%% Commande de la dynamique de rotation
zeta = 0.7;
wn = 20; %[rad/s] !!!
% wn = wn*(180/pi);

Kp_rot = wn^2;
Kd_rot = 2*zeta*wn;

%% Simulation
% Conditions initiales et temps final
z0 = [v_ini, gamma_ini, h_ini, s_ini, theta_ini, q_ini];
tspan = [0, 100];

reltol1 = 10e-10;
options = odeset('reltol', reltol1);
[t, z] = ode45('capsule', tspan, z0, options);
%%
figure
plot(z(:,3), z(:,1))
set(gca,'Xdir','reverse')
title('Vitesse selon altitude de la capsule, sans asservissement')

