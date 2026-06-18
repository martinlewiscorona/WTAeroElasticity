function [Rx, FN, FT, P, a_dyn_out] = BEM_dynamic(v0, omega, pitch, dt, a_dyn_prev)
    % ------------------------------------------------
    % Modified BEM for Task 7: Dynamic Inflow
    % ------------------------------------------------
    B = 3;           
    R = 63;          
    hubrad = 1.5;    
    rou = 1.225;     
    EPS = 0.00001;   
    
    % Wake time constant (approximate inertia of the wake)
    tau = R / v0; 

    BS = importdata('Blade\Blade section\Blade section.dat').data;
    
    % Use persistent variable to cache AD to avoid loading files 3000 times
    persistent AD cached_files
    if isempty(AD) || isempty(cached_files)
        Readfiles = dir(fullfile('Blade\Aero data\','*.dat'));
        cached_files = Readfiles;
        for i=1:length(Readfiles)
            AD{i} = importdata(strcat('Blade\Aero data\',Readfiles(i).name));
        end
    end

    NBS = length(BS);    
    Rx = zeros(NBS, 1); FN = zeros(NBS, 1); FT = zeros(NBS, 1); Mx = zeros(NBS, 1);
    a_dyn_out = zeros(NBS, 1);
    
    for i = 1:NBS
        ADofBS = BS(i,2); 
        r_local = BS(i,3);      
        Rx(i) = r_local;        
        dr = BS(i,4);      
        Theta = BS(i,5);  
        chord = BS(i,6);   
        
        alpha_tbl = AD{ADofBS}(:,1);
        Cl_tbl = AD{ADofBS}(:,2); 
        Cd_tbl = AD{ADofBS}(:,3); 
        Sigma = chord*B / (2*pi*r_local); 
        
        % ITERATIVE SOLVER FOR QUASI-STEADY INDUCTION
        a_qs = 0; a_prime = 0;
        ax = a_qs; ax_prime = a_prime;
        a_qs = ax - 10*EPS; a_prime = ax_prime - 10*EPS; 
        
        numite = 0; 
        while abs(ax - a_qs) >= EPS || abs(ax_prime - a_prime) >= EPS
            numite = numite + 1;
            a_qs = ax;
            a_prime = ax_prime;
            
            Phi_rad = atan((1 - a_qs)*v0 / ((1 + a_prime)*r_local*omega));
            Phi = rad2deg(Phi_rad);
            
            Alpha = Phi - Theta - pitch;
            
            Cla = interp1(alpha_tbl, Cl_tbl, Alpha, 'linear', 'extrap');
            Cda = interp1(alpha_tbl, Cd_tbl, Alpha, 'linear', 'extrap');
            
            Cn = Cla*cosd(Phi) + Cda*sind(Phi);
            Ct = Cla*sind(Phi) - Cda*cosd(Phi);
            
            f_tiploss = B/2*(R - r_local)/(r_local*sind(Phi));
            F_tiploss = (2/pi)*acos(exp(-f_tiploss));
            f_hubloss = B/2*(r_local - hubrad)/(r_local*sind(Phi));
            F_hubloss = (2/pi)*acos(exp(-f_hubloss));
            F = F_tiploss * F_hubloss;
            
            % Guard against division by zero or NaN at tip/root
            if isnan(F) || F < 0.001; F = 0.001; end
            
            ac = 0.2;
            if ax > ac
                K_val = 4*F*sind(Phi)^2 / (Sigma*Cn);
                ax = 0.5*(2 + K_val*(1 - 2*ac) - sqrt((K_val*(1 - 2*ac) + 2)^2 + 4*(K_val*ac^2 - 1)));
            else
                ax = 1 / (4*F*(sind(Phi))^2 / (Sigma*Cn) + 1);
            end
            ax_prime = 1 / (4*F*sind(Phi)*cosd(Phi) / (Sigma*Ct) - 1);
            
            if numite >= 100
                ax = 0.3; ax_prime = 0.1;
                break;
            end
        end
        
        % -------------------------------------------------------------
        % DYNAMIC INFLOW FILTER (FIRST ORDER TIME LAG)
        % Applies inertia to the axial induction factor.
        % -------------------------------------------------------------
        if dt > 0
            % Exponential smoothing filter based on wake time constant
            k_filter = exp(-dt / tau);
            a_dyn = a_dyn_prev(i) * k_filter + ax * (1 - k_filter);
        else
            a_dyn = ax;
        end
        a_dyn_out(i) = a_dyn;
        
        % Recalculate Flow Angle and Forces based on DYNAMIC induction
        Phi_rad_dyn = atan((1 - a_dyn)*v0 / ((1 + ax_prime)*r_local*omega));
        Phi_dyn = rad2deg(Phi_rad_dyn);
        Alpha_dyn = Phi_dyn - Theta - pitch;
        
        Cla_dyn = interp1(alpha_tbl, Cl_tbl, Alpha_dyn, 'linear', 'extrap');
        Cda_dyn = interp1(alpha_tbl, Cd_tbl, Alpha_dyn, 'linear', 'extrap');
        
        Cn_dyn = Cla_dyn*cosd(Phi_dyn) + Cda_dyn*sind(Phi_dyn);
        Ct_dyn = Cla_dyn*sind(Phi_dyn) - Cda_dyn*cosd(Phi_dyn);
        
        % Calculate unsteady aerodynamic forces
        V_rel_sq = (r_local*omega*(1 + ax_prime))^2 + (v0*(1 - a_dyn))^2;
        FN(i) = 0.5 * rou * V_rel_sq * chord * Cn_dyn * dr;
        FT(i) = 0.5 * rou * V_rel_sq * chord * Ct_dyn * dr;
        Mx(i) = FT(i) * r_local;
    end
    
    M_torque = sum(Mx); 
    P = M_torque * omega * B * 0.944; 
end