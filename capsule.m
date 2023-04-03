function f = capsule(t, z)
    ctl = 0; % 0: pas d'asservissement, 1: asservissement
    load('data.mat')
    % Variables d'état
    v = z(1);
    gamma = z(2);
    h = z(3);
    s = z(4);
    theta = z(5);
    q = z(6);
    
    rho = rho0*exp(-h/hs);
    r = (R_mars+h);
    r_fin = R_mars + h_fin;
    g_r = mu_mars/r^2;
    Pdyn = 0.5*rho*v^2; 
    
    alpha = theta - gamma;
    
    if ctl == 0
        delta_cmd = 0.0;
    else
        % Partie asservissement translation
        delta_v_aero = v_fin1 - sqrt(v^2 + 2*mu_mars*(1/r_fin - 1/r));
        gamma_ref = asind((0.5 * B *hs *(rho_fin- rho))/(log(1+delta_v_aero/v)));
        g_theta = (Pdyn*S*CLa)/(v*m);
        theta_eq = -gamma + ((cosd(gamma)*v*m)/(Pdyn*S*CLa))*(v^2/(R_mars+h)^2 - g_r);
        theta_cmd = theta_eq  + (Kp_tra/g_theta)*(gamma_ref - gamma);
        
        alpha = theta_cmd - gamma;
        
        % S'assurer que theta cmd ne dépasse pas 60 deg
        % aka ajouter bout de code ici
        
        % Partie asservissement rotation
        delta_eq = -((CMa*(theta-gamma) + (d/(2*v))*CMq*q))/CMd;
        g_delta = Pdyn*S*d*CMd/J;
        delta_cmd = delta_eq + Kp_rot*(theta_cmd - theta)/g_delta + Kd_rot*q/g_delta;
         
    end
    
    % Parties utilisés pour les équations d'ordre 1 finales
    Daero = Pdyn * S * CD0;
    Maero = Pdyn*S*d*(CMa*alpha + d*CMq*q/(2*v) + CMd*delta_cmd);
    Laero = Pdyn*S*CLa*alpha;
    
    % v dot
    f(1) = -Daero/m - g_r*sind(gamma);
    
    % gamma dot
    f(2) = (1/v)*(Laero/m + (v^2/r - g_r)*cosd(gamma));
    
    % h dot
    f(3) = v * sind(gamma);
    
    % s dot
    f(4) = (v/r) * cosd(gamma);
    
    % theta dot 
    f(5) = q;
    
    % q dot
    f(6) = (1/J) * Maero;
    
    % Daero > 2650
%     f(7) = Daero > D_aero_max;
    
    f = f';


